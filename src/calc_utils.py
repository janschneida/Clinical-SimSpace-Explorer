import pandas as pd
import numpy as np
import sys
import os
import functools
from pygosemsim import term_set
from pygosemsim import similarity
import multiprocessing
import concurrent.futures
from scipy.optimize import linear_sum_assignment
sys.path.append(os.getcwd())
import src.data_utils as du
from src.taxodist import td_calc 
import math
import time
from scipy.spatial.distance import euclidean
import dask
import dask.distributed as dd
import csv

import time

dist_scores = []
mut_prof_dists = []
dem_dists = []
icds_dists = []
mut_prof_sims_lookup = {}
mut_prof_sims_filtered_lookup = {}

def getDistScores(df: pd.DataFrame,weights:dict = {},calc_dem_dist = False,calc_icds_dist=False,calc_mut_dist = False,term_source='panther',filter_genes=False) -> pd.DataFrame:
    ''' magical method that computes similarity scores
        returns normalized, weighted & merged pairwise distance matrix
    '''   
    global mut_prof_dists
    global dem_dists
    global icds_dists
    global dist_scores
    # check if we need to recalculate the dems (is the case if a dem feature-weight changed or on init)
    if calc_dem_dist:
        dem_start = time.time()
        dem_dists = __getDemographicDist__(df,weights)
        dem_end = time.time()
        dem_time = dem_end - dem_start
        print('finished dem dists',dem_dists.shape,dem_time)
    
    if calc_icds_dist:
        icd_start = time.time()
        icds_dists = __getICDsDist__(df)
        icd_end = time.time()
        icd_time = icd_end - icd_start
        print('finished icd dists',icds_dists.shape,icd_time) 
        
    # check if we need to (re)calculate the whole matrix or just change the weight
    if calc_mut_dist:
        mut_start = time.time()
        mut_prof_dists = 1 - __getMutationProfileSims__(du.getMutationProfileDictFromDataFrame(df),term_source=term_source,filter_genes=filter_genes)
        mut_end = time.time()
        mut_time=mut_end-mut_start
        print('finished mut dists',mut_prof_dists.shape,mut_time)

    # NOTE hook for SNF
    dist_scores = ((dem_dists + weights['mutations']*mut_prof_dists) + weights['icds']*icds_dists )/3
    # NOTE hook to change MDS to t-SNE
    df[['x','y']] = du.getMDSMatrix(dist_scores)
    # df[['x','y']] = du.getTSNEMatrix(dist_scores)
    return df

def __getICDsDist__(df: pd.DataFrame) -> np.ndarray:
    td = td_calc.Taxodist()
    icd_sets = list(du.getICDsDictFromDataFrame(df).values())
    # TODO which icd version?
    tree = du._initICDTree_()
    matrix = td.calc_set_sim(
        sets=icd_sets,
        tree=tree,
        ic_mode='levels',cs_mode='leacock_chodorow',setsim_mode='cosine',
        #ic_mode='levels',cs_mode='leacock_chodorow',setsim_mode='bipartite_matching',
        normalize=True)
    return 1 - matrix

def __getDemographicDist__(df: pd.DataFrame, feature_weights: dict  = {}) -> np.ndarray:
    '''
        method that calculates euclidian distance for numerical features
        returns normalized pairwise distance matrix
    '''
    
    demographics = du.__getDemographicsDataFrame__(df)
    dist_matrix = np.zeros((len(df), len(df)))
    
    for i in demographics.index:
        pat1_vector = list(demographics.iloc[i])
        for j in demographics.index:
            pat2_vector = list(demographics.iloc[j])
            # Calculate the weighted Euclidean distance
            dist = 0.0
            for index, value1 in enumerate(pat1_vector):
                col = demographics.columns[index]
                weight = feature_weights.get(col, 1.0)  # Default weight is 1.0 if not specified
                value2 = pat2_vector[index]
                if not np.isnan(value1) and not np.isnan(value2):
                    dist += weight * (value1 - value2) ** 2
                dist = np.sqrt(dist)
            
            dist_matrix[i,j] = dist
            dist_matrix[j,i] = dist_matrix[i,j]
            
    if np.max(dist_matrix) == 0:
        return dist_matrix
    else:
        return (dist_matrix-np.min(dist_matrix))/(np.max(dist_matrix)-np.min(dist_matrix))

def __getMutationProfileSims__(patientMutations: dict,term_source='panther',scale_to_muts=True,scale_to_go_terms=False,filter_genes=False) -> np.ndarray:
    '''
    Returns pairwise-patient-similarity-matrix based on mutation profiles.
    patientMutations is a dict that maps patientId's to lists of genes
    '''
    
    global mut_prof_sims_lookup
    global mut_prof_sims_filtered_lookup
    if filter_genes:
        mut_prof_sims = mut_prof_sims_filtered_lookup.copy()
    else:
        mut_prof_sims = mut_prof_sims_lookup.copy()
    du._initPyGoGO_()
    term_dict = {}
    uniques = du.getUniqueMutations(patientMutations)
    # TODO from files!
    if term_source == 'panther':
        term_dict = du.__initPantherDBtermDict__(uniques)
    elif term_source == 'pygosemsim':
        term_dict = du.__initPyGOSemSimTermDict__(uniques)

    # precompute called every time in case new genes are included (simply adds those to the lookup)
    # NOTE those files would be different depending on the algorithm
    __initGeneSimLookUp__(uniques,term_dict,scale_to_go_terms)
    keysList = list(patientMutations.keys())
    matrix = np.zeros(shape=(len(patientMutations),len(patientMutations)))
    for i, patMuts1 in enumerate(keysList):
        for j, patMuts2 in enumerate(keysList[i:]):
            # check if mutSim has already been calculated before
            if not (patMuts1,patMuts2) in mut_prof_sims and not (patMuts2,patMuts1) in mut_prof_sims:
                mutSim = __getMutProfileSim__(patientMutations[patMuts1],patientMutations[patMuts2],scale_to_muts)
                # TODO wieder reinnehmen
                mut_prof_sims[(patMuts1,patMuts2)] = mutSim
                mut_prof_sims[(patMuts2,patMuts1)] = mutSim
            matrix[i, j+i] = mut_prof_sims.get((patMuts1,patMuts2))
            matrix[j+i, i] = mut_prof_sims.get((patMuts1,patMuts2))
    # update lookups
    if filter_genes:
        mut_prof_sims_filtered_lookup.update(mut_prof_sims)
    else:
        mut_prof_sims_lookup.update(mut_prof_sims)
    # scale matrix
    return (matrix-np.min(matrix))/(np.max(matrix)-np.min(matrix))

def __getMutProfileSim__(muts1: list,muts2: list,scale_to_muts=True):
    '''
    Calculates similarity of two mutation profiles. 
    This is done by calculalating the set similarity of the respective gene-lists.
    For this we need to calculate pairwise gene-similarities.
    The gene similarities on the other hand consist of GO-term-set similarities.
    
    We currently use the bipartite graph matching algorithm for the mut-profile set-sim without any scaling 
    (might not be needed..some analysis of the mutation set size could be useful tho)
    
    yikes brotheeeer thats a lot of calculating right here :'))))))
    '''
    gene_sim_matrix = np.zeros(shape=(len(muts1),len(muts2)))

    # NOTE removing duplicates
    # NOTE we might also use indexed search but lets keep it like this for now
    muts1 = list(set(muts1))
    muts2 = list(set(muts2))
    
    gene1_index = 0
    for gene1 in muts1:
        gene1_index = muts1.index(gene1)
        gene2_index = 0    
        for gene2 in muts2:
            gene2_index = muts2.index(gene2)
            # calculate similarity of two genes (function does GO term scaling)
            geneSim = du.gene_sim_look_up.get((gene1,gene2))
            # TODO wo kommen die Nones her???
            
            # this little line saves us lots of space
            if geneSim is None:
                geneSim = 0.0
            gene_sim_matrix[gene1_index,gene2_index] = geneSim

    # TODO add other set sim possibilities
    # TODO account for empty 
    # use kuhn munkres algorithm 
    row_ind, col_ind = linear_sum_assignment(cost_matrix=gene_sim_matrix,maximize=True)
    # calculate mut profile sim
    mutSim = gene_sim_matrix[row_ind, col_ind].sum() 
    if scale_to_muts:
        mutSim = __getScaledSetSim__(mutSim,len(muts1),len(muts2))
        
    return mutSim

def __getGeneSim__(genePair,terms_dict:dict,go, scale = False):
    """ Returns scaled GO-term-set-sim for the given genes. """
    gene1, gene2 = genePair
    terms1 = terms_dict.get(gene1)
    terms2 = terms_dict.get(gene2)
    # NOTE adjust return values or always scale - depending on algorithm
    if not terms1 or not terms2:
        return gene1, gene2, 0.0
    # if terms1 == terms2:
    #     return gene1, gene2, 1.0
    
    # NOTE for now we keep Resnik + BMA 
    sf = functools.partial(term_set.sim_func,go, similarity.resnik)
    gene_sim = term_set.sim_bma(terms1, terms2, sf)
    if scale:
        gene_sim = __getScaledSetSim__(gene_sim,len(terms1),len(terms2))
    return gene1,gene2,gene_sim

def __getScaledSetSim__(setSim, len1, len2):
    ''' Used to scale the set-similarities to account for differences in set-sizes that might impair the accuracy of the calculations. '''
    setDiff = abs(len1 - len2)
    maxSim = min(len1,len2)
    return setSim/(maxSim + math.log(setDiff + 1))

def __initGeneSimLookUp__(genes: list,term_dict: dict, scale_to_go_terms=False,recalc=False):
    for gene1 in genes:
        for gene2 in genes:
            if (gene1, gene2) not in du.gene_sim_look_up and (gene2, gene1) not in du.gene_sim_look_up:
                #genePair,peneJair,gene_sim = __getGeneSim__((gene1, gene2), term_source, scale_to_go_terms)
                gene1,gene2,gene_sim = __getGeneSim__((gene1, gene2),term_dict, du.go, scale_to_go_terms)
                # dont save 0.0s because we handle 
                if gene_sim == 0.0:
                    continue
                du.gene_sim_look_up[(gene1,gene2)] = gene_sim
                du.gene_sim_look_up[(gene2,gene1)] = gene_sim
    
# NOTE this does not work with dash warum auch immer
# def __parProc_precomputeGeneSimilarities__(genes: list, term_dict:dict,go, scale_to_go_terms=False):
#     gene_sim_look_up = {}

#     # Generate all gene pairs to be processed
#     gene_pairs = [(gene1, gene2) for gene1 in genes for gene2 in genes if (gene1, gene2) not in gene_sim_look_up and (gene2, gene1) not in gene_sim_look_up]

#     # Create a process pool with a suitable number of processes
#     with multiprocessing.Pool() as pool:
#         # Map the tasks to the process pool for parallel execution
        
#         results = pool.starmap(__getGeneSim__, [(pair,term_dict,go,scale_to_go_terms) for pair in gene_pairs])

#     for result in results:
#         print(result)
#     sim_results = []
#     with concurrent.futures.ProcessPoolExecutor() as executor:
#         # Stage 1: Compute similarities within each batch
#         for gene_pair in gene_pairs:
#             sim_results.append(executor.submit(__getGeneSim__, gene_pair, term_source, scale_to_go_terms))
    
#         for future in concurrent.futures.as_completed(sim_results):
#             print(future.result())

#     du.gene_sim_look_up.update(gene_sim_look_up)

# # NOTE batches also dont offer any speedup :')
# def __parThreads_batches_precomputeGeneSimilarities__(genes: list, term_source='panther', scale_to_go_terms=False, batch_size=100):
#     gene_sim_look_up = {}

#     # Split genes into batches
#     gene_batches = [genes[i:i+batch_size] for i in range(0, len(genes), batch_size)]

#     # Create a ThreadPoolExecutor with a suitable number of workers
#     with concurrent.futures.ThreadPoolExecutor() as executor:
#         # Stage 1: Compute similarities within each batch
#         batch_results = []
#         for gene_batch in gene_batches:
#             batch_results.append(executor.submit(calculate_batch_similarity, gene_batch, term_source, scale_to_go_terms))

#         # Collect the results of Stage 1
#         for future in concurrent.futures.as_completed(batch_results):
#             batch_similarities = future.result()
#             for gene_pair, gene_sim in batch_similarities:
#                 gene_sim_look_up[gene_pair] = gene_sim

#         # Stage 2: Compute similarities across different batches
#         cross_batch_results = []
#         for i in range(len(gene_batches)):
#             for j in range(i+1, len(gene_batches)):
#                 batch1 = gene_batches[i]
#                 batch2 = gene_batches[j]
#                 cross_batch_results.append(executor.submit(calculate_cross_batch_similarity, batch1, batch2, term_source, scale_to_go_terms))

#         # Collect the results of Stage 2
#         for future in concurrent.futures.as_completed(cross_batch_results):
#             cross_batch_similarities = future.result()
#             for gene_pair, gene_sim in cross_batch_similarities:
#                 gene_sim_look_up[gene_pair] = gene_sim

#     du.gene_sim_look_up.update(gene_sim_look_up)
#     #return gene_sim_look_up

# def calculate_batch_similarity(genes_batch, term_source, scale_to_go_terms):
#     batch_similarities = []

#     for gene1 in genes_batch:
#         for gene2 in genes_batch:
#             if gene1 != gene2:
#                 gene_sim = __getGeneSim__(gene1, gene2, term_source, scale_to_go_terms)
#                 batch_similarities.append(((gene1, gene2), gene_sim))

#     return batch_similarities

# def calculate_cross_batch_similarity(batch1, batch2, term_source, scale_to_go_terms):
#     cross_batch_similarities = []

#     for gene1 in batch1:
#         for gene2 in batch2:
#             gene_sim = __getGeneSim__(gene1, gene2, term_source, scale_to_go_terms)
#             cross_batch_similarities.append(((gene1, gene2), gene_sim))

#     return cross_batch_similarities


# # NOTE the threads dont offer any speedup (??) so very much useless
# def __parThreads_noBatches_precomputeGeneSimilarities__(genes: list, term_source='panther', scale_to_go_terms=False):
#     gene_sim_look_up = {}

#     # Generate all gene pairs to be processed
#     gene_pairs = [(gene1, gene2) for gene1 in genes for gene2 in genes if (gene1, gene2) not in gene_sim_look_up and (gene2, gene1) not in gene_sim_look_up]

#     # Create a ThreadPoolExecutor with a suitable number of workers
#     with concurrent.futures.ThreadPoolExecutor() as executor:
#         # Submit tasks to the thread pool for parallel execution
#         results = [executor.submit(__getGeneSim__, pair, term_source, scale_to_go_terms) for pair in gene_pairs]

#         # Collect the results as they complete
#         for future in concurrent.futures.as_completed(results):
#             result_pair = future.result()
#             gene_sim_look_up[result_pair[0][0]] = result_pair[0][1]
#             gene_sim_look_up[result_pair[1][0]] = result_pair[1][1]

#     du.gene_sim_look_up.update(gene_sim_look_up)

# # NOTE maybe thissssssssss???
# import dask
# import dask.bag as db
# from dask.distributed import Client

# def __Dask_precomputeGeneSimilarities__(genes: list, terms_dict,go, scale_to_go_terms=False):
#     gene_sim_look_up = {}

#     # Convert the gene list to a Dask bag
#     genes_bag = db.from_sequence(genes)
#     gene_pairs = [(gene1, gene2) for gene1 in genes for gene2 in genes if (gene1, gene2) not in gene_sim_look_up and (gene2, gene1) not in gene_sim_look_up]
#     gene_pairs_bag = db.from_sequence(gene_pairs,partition_size=1)

#     # Compute the gene similarities in parallel using Dask
#     # gene_similarities = genes_bag.product(genes_bag).map(
#     #     lambda x: __getGeneSim__(x, getGenePairTerms(x,term_source), scale_to_go_terms)).compute()
    
#     #if term_source == 'panther':
#     gene_similarities = gene_pairs_bag.map(
#         lambda x: __getGeneSim__(x,terms_dict,go,scale_to_go_terms)).compute()
#     #if term_source == 'pygosemsim':
#         # gene_similarities = genes_bag.product(genes_bag).map(
#         # lambda x: __getGeneSim__(x,terms_dict, scale_to_go_terms)).compute()


#     # Collect the results and update the gene_sim_look_up dictionary
#     for result in gene_similarities:
#         gene1, gene2, gene_sim = result
#         gene_sim_look_up[(gene1,gene2)] = gene_sim
#         gene_sim_look_up[(gene2,gene1)] = gene_sim

#     du.gene_sim_look_up.update(gene_sim_look_up)
    
# def getGenePairTerms(genePair:tuple,term_source):
#     gene1, gene2 = genePair
#     try:
#         if term_source == 'panther':
#             terms1 = du.panther_terms.get(gene1)
#             terms2 = du.panther_terms.get(gene2)
#         elif term_source == 'pygosemsim':
#             terms1 = du.pygosemsim_terms.get(gene1)
#             terms2 = du.pygosemsim_terms.get(gene2)
#         else:
#             raise ValueError("Unsupported term_source: ", term_source)
#     except ValueError as err:
#         print(err.args)
#         sys.exit()
#     return (terms1,terms2)

# def fillLookUps(pat_mutations: dict):
#     '''
#     precomputes similarity scores of all unique mutations in the cohort
#     '''
#     # TODO might be usefull to move to the start of calcMutProfile() 
#     # instead of calling the getUnique and init
#     unique_muts = du.getUniqueMutations(pat_mutations)
#     # with concurrent.futures.ProcessPoolExecutor() as executor:
#     #     for i, patMuts1 in enumerate(keysList):
#     #         # Submit each pairwise calculation as a separate task
#     #         results = [executor.submit(calcPairwiseMutProfileSim, patientMutations[patMuts1], patientMutations[patMuts2], scale_to_muts, scale_to_go_terms)
#     #                     for j, patMuts2 in enumerate(keysList[i:])]

#     return

if __name__ == '__main__':
    multiprocessing.freeze_support()