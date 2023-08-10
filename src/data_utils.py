import json
import requests
import sys
import os
from pygosemsim import download
from pygosemsim import similarity
from pygosemsim import graph
from cbio_py import cbio_mod as cb
# is problem, "needs" network connection
import pandas as pd
sys.path.append(os.getcwd())
import src.gaf_annotation as gaf_annotation
from unipressed import IdMappingClient
import re
import time
import numpy as np
from sklearn.manifold import MDS, TSNE
from multiprocessing import freeze_support
import csv
from src.taxodist import tree_parsers

import time

# globals
annots = []
go = []
gene_sim_look_up = {}
pygosemsim_terms = {}
panther_terms = {}
uniprot_lookup = {}
oncoKB_lookup = []
tree = []

def getGOAnnotationsFromPantherDB(genes: list, ontology:str='BP') -> dict:
    '''
    Returns dict of GO annotations for a given set of genes from given ontology (BP,CC,MF)
    using the pantherDB API
    
    "Each identifier to be delimited by comma i.e. ',. 
    Maximum of 1000 Identifiers can be any of the following: 
    Ensemble gene identifier, Ensemble protein identifier, 
    Ensemble transcript identifier, Entrez gene id, gene symbol, 
    NCBI GI, HGNC Id, International protein index id, NCBI UniGene id, UniProt accession andUniProt id"
    
    possible ontologies: BP, CC, MF
    
    '''
    genes_string = ','.join(genes)
    url = "http://pantherdb.org/services/oai/pantherdb/geneinfo?geneInputList={genes}&organism=9606".format(genes=genes_string)
    response = requests.post(url)
    data = json.loads(response.text)

    # NOTE using GO SLIMs natively! 
    content = 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_{ontology}'.format(ontology=ontology.upper())

    # retrieve data from data
    gene_dict = {}
    # iterate over genes
    for i, gene in enumerate(data['search']['mapped_genes']['gene']):
        go_terms = []
        # TODO crashes when only one gene is queried but ehehh lets do that later
        if gene.get('annotation_type_list','nope') == 'nope':
            continue
        annot_lists = gene['annotation_type_list']['annotation_data_type']
        for annot_list in annot_lists:
            if type(annot_list) == str:
                break
            # find desired annotation list in data
            if annot_list['content'] == content:
                # append annotated terms
                for annot in annot_list['annotation_list']['annotation']:
                    # if only 1 term is annotated (weird workaround) but eh
                    if type(annot) is str:
                        id = annot_list['annotation_list']['annotation']['id']
                        go_terms.append(id)
                        break
                    else:
                        go_terms.append(annot['id'])
        # map genes from list to terms
        # NOTE we dont need to save empty lists, when we use .get() on the gene_dict
        # NOTE is this better? or maybe save empty list to ensure we dont query genes we havent seen before
        if go_terms:
            uniprotID = __getAccession__(gene)
            gene_name = uniprot_lookup.get(uniprotID)
            if gene_name:
                gene_dict[gene_name] = go_terms
    
    return gene_dict

def filterOncoGenesInPatientDataFrame(df:pd.DataFrame):
    filtered_df = df.copy()
    filtered_df['mutations'] = filtered_df['mutations'].apply(_filter_mut_prof_)
    return filtered_df

def _filter_mut_prof_(cell_value: str):
    values = cell_value.split()  # Split the cell value into a list of values
    filtered_values = [value for value in values if value in oncoKB_lookup]  # Filter values based on the lookup list
    return ' '.join(filtered_values)  # Join the filtered values back into a string  

def __getAccession__(gene:dict) -> str:
    accession = gene['accession']
    result = re.search(r"UniProtKB=(.*)", accession)
    return result.group(1)

def getGOAnnotationsFromPyGoSemSim(gene, aspect='P') -> list:
    '''
    Returns list of GO terms from annots based on the HUGO gene Symbol
    the aspect variable can be used to filter the annotations based on their corresponding ontology: \n
    P (biological process) F (molecular function) C (cellular component)
    '''
    terms = []
    for gene_iter in annots:
        gene_iter = annots[gene_iter]
        if gene_iter['db_object_symbol'] == gene:
            annotations = gene_iter['annotation']
            for term in annotations:
                if annotations[term]['aspect'] == aspect:
                    terms.append(term)
    return terms

def getMutationProfileDictFromStudy(studyId, groupByFeature='gene',geneFeature='hugo') -> dict:
    '''
    returns dict that maps patientIDs out of cBioPortal study with 'studyId' to mutation profiles
    - groupByFeature allows to retrieve the mut_dict to be grouped 
      by the patientId and mapped to a certain attr from muts
    '''
    mut_dict = {}
    # NOTE append == 'yes' as default adds _all & _mutations to the studyId for the molecularProfileId & sampleListId
    # NOTE can be changed by setting 'append' to literally anything else if we would want to access other stuff
    muts = cb.getMutationsInMolecularProfile(molecularProfileId=studyId,sampleListId=studyId)

    for mutation in muts:
        patient_id = mutation['patientId']
        # TODO maybe allow other attr but right now lets keep it a little simple and stick to HUGO symbols
        if groupByFeature == 'gene':
            gene = mutation[groupByFeature].hugoGeneSymbol
        else :
            gene = mutation[groupByFeature]
        # start gene list by grouping by patientId or append genes to gene list
        mut_dict.setdefault(patient_id, []).append(gene)
    
    return mut_dict

def minMaxScaling(df: pd.DataFrame):
    return (df - df.min()) / (df.max() - df.min())

def __initPyGoAnnots__():
    ''' Inits the annots look up for the PyGOSemSim GO-terms'''
    global annots 
    if annots:
        return
    # annots = annotation.from_resource('goa_human')
    annots = gaf_annotation.from_resource('goa_human')

def _initPyGoGO_():
    ''' Inits the GO object needed for sim calculations.'''
    # TODO could also use the most up-to-date obo every time:
    # go_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
    # urllib.request.urlretrieve(go_url, 'go-basic.obo')
    global go 
    if go:
        return
    go = graph.from_obo('resources/application/go-basic.obo')
    similarity.precalc_lower_bounds(go)
 
def __initPantherDBtermDict__(genes: list, ontology = 'BP',reload=False):
    global panther_terms
    start = time.time()
    if not reload:
        with open('resources/application/panther_term_dict.csv', 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                gene = row['Gene']
                go_terms = row['GO Terms'].split(' ')
                panther_terms[gene] = go_terms
    
    if panther_terms:
        return panther_terms
    __initUniprotLookUp__(genes)
    panther_terms = {}
    if len(genes) < 1000:
        panther_terms = getGOAnnotationsFromPantherDB(genes, ontology)
    else:
        batch_size = 999  # Set the batch size to 999 to ensure it is less than 1000
        num_batches = len(genes) // batch_size

        for i in range(num_batches):
            batch_genes = genes[i * batch_size: (i + 1) * batch_size]
            batch_terms = getGOAnnotationsFromPantherDB(batch_genes, ontology)
            panther_terms.update(batch_terms)

        remaining_genes = genes[num_batches * batch_size:]
        if remaining_genes:
            remaining_terms = getGOAnnotationsFromPantherDB(remaining_genes, ontology)
            panther_terms.update(remaining_terms)
    end = time.time()
    gene_load_time = end-start
    print('time to load',len(genes),'genes from panther:',gene_load_time)
    return  panther_terms

def _initOncoKBLookUp_(reload=False):
    global oncoKB_lookup
    
    # no reload -> load from file
    if not reload:
        with open('resources/application/oncoKB_lookup.csv', 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                gene = row['Oncogene/TSG Gene']
                oncoKB_lookup.append(gene)

    url = 'https://www.oncokb.org/api/v1/utils/allCuratedGenes?includeEvidence=true'
    response = requests.get(url)
    data = json.loads(response.text)
    for gene in data:
        if gene['oncogene'] or gene['tsg']:
            oncoKB_lookup.append(gene['hugoSymbol'])

def __initPyGOSemSimTermDict__(genes: list, ontology ='P',reload=False):
    global pygosemsim_terms
    
    if not reload:
        with open('resources/application/pygosemsim_term_dict.csv', 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                gene = row['Gene']
                go_terms = row['GO Terms'].split(' ')
                pygosemsim_terms[gene] = go_terms
    
    if pygosemsim_terms:
        return pygosemsim_terms
    __initPyGoAnnots__()
    for gene in genes:
        pygosemsim_terms[gene] = getGOAnnotationsFromPyGoSemSim(gene,ontology)
    return pygosemsim_terms

def getFilteredWeights(df: pd.DataFrame, feature_weights: dict) -> dict:
    '''
    Sets feature weights to zero for not-included features
    '''
    filteredWeights = {}
    for feature, weight in feature_weights.items():
        break
    return filteredWeights

def getUniqueMutations(pat_mutations: dict) -> list:
    pat_mutations = list(pat_mutations.values())
    pat_mutations_flat = list()

    for sub_list in pat_mutations:
        pat_mutations_flat += sub_list
        
    return list(set(pat_mutations_flat))

def __initUniprotLookUp__(gene_names: list,reload=False):
    '''inits the uniprotID to gene-name look up needed to map queried genes to their terms'''
    
    global uniprot_lookup
    # no reload -> load from file
    if not reload:
        with open('resources/application/uniprot_lookup.csv', 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                gene = row['Uniprot']
                hugo_name = row['HUGO Gene Name']
                uniprot_lookup[gene] = hugo_name
        return

    request = IdMappingClient.submit(
        source="GeneCards", dest="UniProtKB", ids=gene_names
    )
    time.sleep(10)
    for res in list(request.each_result()):
        uniprot_lookup[res['to']]=res['from']
    print('finished uniprot lookup')
        
def flushGeneSimLookUp():
    '''
    Method to call when term_source parameter changes
    '''
    global gene_sim_look_up
    gene_sim_look_up = {}

def getMutationProfileDictFromDataFrame(df: pd.DataFrame) -> dict:
    patdict = {}
    for values in df[['patient_id','mutations']].values:
        pat_id = values[0]
        mut_profile = values[1]
        mut_profile = str(mut_profile).split(' ')
        patdict[pat_id]=mut_profile
    return patdict

def getICDsDictFromDataFrame(df: pd.DataFrame) -> dict:
    patdict = {}
    for values in df[['patient_id','icds']].values:
        pat_id = values[0]
        icd_codes = values[1]
        icd_codes = str(icd_codes).split(' ')
        patdict[pat_id]=icd_codes
    return patdict

def getAllStudyIDsFromCBioPortal():
    studies = cb.getAllStudies()
    study_ids = []
    for study in studies:
        study_ids.append(study['studyId'])
    return study_ids

def __getDemographicsDataFrame__(df: pd.DataFrame) -> pd.DataFrame:
    '''
    returns normalized demographics DataFrame
    '''
    df_features = df.columns
    demographics = pd.DataFrame()
    if 'age' in df_features:
        df['age'] = df['age'].apply(pd.to_numeric)
        min_age = df['age'].dropna().min()
        max_age = df['age'].dropna().max()

        # Normalize the values in the 'age' column between 0 and 1
        demographics['age'] = df['age'].map(lambda x: (x - min_age) / (max_age - min_age) if pd.notnull(x) else x)

    if 'age at diagnosis' in df_features:
        df['age at diagnosis'] = df['age at diagnosis'].apply(pd.to_numeric)
        min_age_at_diag = df['age at diagnosis'].dropna().min()
        max_age_at_diag = df['age at diagnosis'].dropna().max()

        demographics['age at diagnosis'] = df['age at diagnosis'].map(lambda x: (x - min_age_at_diag) / (max_age_at_diag - min_age_at_diag) if pd.notnull(x) else x)

    if 'sex' in df_features:
        demographics['sex'] = df['sex']
        demographics = demographics.replace(['M','Male'],0)
        demographics = demographics.replace(['F','Female'],1)
    return demographics


def getMDSMatrix(dist_matrix: np.ndarray) -> pd.DataFrame:
    """Computes multi-dimensionally-scaled two-dimensional concept-coordinates based on a pairwise-distance-matrix"""
    # use MDS to compute the relative distances of the distinct concepts
    embedding = MDS(n_components=2,dissimilarity='precomputed')
    return pd.DataFrame(embedding.fit_transform(dist_matrix),columns=['x','y'])

def getTSNEMatrix(dist_matrix: np.ndarray) -> pd.DataFrame:
    """Computes multi-dimensionally-scaled two-dimensional concept-coordinates based on a pairwise-distance-matrix"""
    # use MDS to compute the relative distances of the distinct concepts
    embedding = TSNE(n_components=2)
    return pd.DataFrame(embedding.fit_transform(dist_matrix),columns=['x','y'])

def getDataFrameFromCBioPortalStudies(studies_list:list) -> pd.DataFrame:
    '''
    Returns DataFrame that holds combined information of all patients of given studies.\n
    Currently only includes features all studies have in common and drops the others. 
    '''
    dfs = []
    for study in studies_list:
        dfs.append(_getDataFrameFromCBioPortalStudy_(study))
    
    df = pd.concat(dfs, ignore_index=True)
    # NOTE we are currently only including features that all studies share
    # that means, we drop the NaN containing columns
    return df#.dropna(axis=1)

def _initICDTree_():
    global tree
    if not tree:
        tree = tree_parsers.getICD10GMTree()
    return tree

def _getDataFrameFromCBioPortalStudy_(studyId:str):
    # check if study is saved locally
    start = time.time()
    read_dir = 'resources/local_studies'
    local_studies = os.listdir(read_dir)
    if studyId in local_studies:
        df = pd.read_csv(read_dir+'/'+studyId)
        end = time.time()
        load_public_study_time = end-start
        print('time to load public study:', studyId,':',load_public_study_time)
        return df.drop(labels='Unnamed: 0',axis=1)
    
    attributes = ['SEX','AGE','AGE_AT_DIAGNOSIS','ICD_10']
    attr_dict = {'SEX':'sex','AGE':'age','AGE_AT_DIAGNOSIS':'age at diagnosis','ICD_10':'icds'}
    pat_list = []
    pats = cb.getAllPatientsInStudy(studyId)
    for pat in pats:
        patient_id = pat['patientId']
        row = {}
        row['cohort'] = studyId
        row['patient_id'] = patient_id
        pat_data = cb.getAllClinicalDataOfPatientInStudy(studyId=studyId,patientId=patient_id)
        for attr in pat_data:
            if attr['clinicalAttributeId'] in attributes:
                row[attr_dict[attr['clinicalAttributeId']]] = attr['value']
        pat_list.append(row)
        
    # add mutations DataFrame
    muts_dict = {}
    for patient, genes in getMutationProfileDictFromStudy(studyId).items():
        muts_dict[patient] = ' '.join(genes)
    muts_df = pd.DataFrame(list(muts_dict.items()), columns=['patient_id', 'mutations'])
    end = time.time()
    load_public_study_time = end-start
    print('time to load public study:', studyId,':',load_public_study_time)
    return pd.merge(pd.DataFrame(pat_list),muts_df, on='patient_id')

def loadGeneSimLookUp(term_source='panther'):
    ''' Loads precalculated gene similarities into geneSimLookUp to enhance mutSim calculation.
        Can be extended to load as many pre calculated gene similarities as desired
    '''
    global gene_sim_look_up
    if term_source == 'panther':
        with open('resources/application/oncokb_gene_sim_lookup.csv', 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                genePair = row['Gene Pair'].split(' ')
                geneSim = row['Gene Sim']
                if geneSim:
                    geneSim = float(geneSim)
                else:
                    # TODO sort out when to save the 0.0s/when not to and if that might be a problem
                    geneSim = 0.0
                gene_sim_look_up[(genePair[0],genePair[1])] = geneSim
                gene_sim_look_up[(genePair[1],genePair[0])] = geneSim
                # TODO we do not have to save 0.0 values - we can just impute them when reading and i guess were a lread ydoing so


if __name__ == '__main__':
    freeze_support()