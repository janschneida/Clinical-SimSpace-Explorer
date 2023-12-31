import sys
import os
import math
import random
import xml.etree.ElementTree as ET
sys.path.append(os.getcwd())

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import treelib
from pandas.core.frame import DataFrame
from scipy.spatial import distance_matrix
from sklearn.manifold import MDS
from treelib.node import Node
from treelib.tree import Tree
from src.taxodist import cs_algorithms
from src.taxodist import ic_algorithms
from src.taxodist import setsim_algorithms
from numpy import ndarray

# global hashtables
ic_table = {}
cs_table = {}
depth = 0

def iterateOverDiags(parent: ET.Element, parent_node: Node, tree: Tree):
    for diag in parent.iter('diag'):
        diag_name = diag.find('name').text
        if not tree.contains(diag_name):
            diag_node = tree.create_node(diag_name, diag_name, parent=parent_node)
            iterateOverDiags(diag,diag_node,tree)

def getIC(concept: str, tree: Tree, ic_mode: str):
    """
    Returns information content of a given concept 
    based on the IC algorthms from https://doi.org/10.1186/s12911-019-0807-y
    
    """
    global ic_table 
    
    if ic_table.get(concept):
        return ic_table[concept]
    
    try:
        if ic_mode == 'levels':
            # IC calculation based on Boriah et al. https://doi.org/10.1137/1.9781611972788.22 
            ic = tree.depth(concept)
        elif ic_mode == 'sanchez':
            ic = ic_algorithms.getICSanchez(concept,tree)
        else:
            raise ValueError('Unsupported IC-mode: ',ic_mode)
        ic_table[concept] = ic
        return ic
    except ValueError as err:
        print(err.args)
        sys.exit()

def getLCA(concept1: str, concept2: str, tree: Tree, ic_mode: str) -> str:
    """Return lowest common ancester of two concepts."""
    lca = 0
    ca = list(getAncestors(concept1, tree).intersection(getAncestors(concept2, tree)))
    if len(ca) != 0:
        lca = ca[0]
        for concept in ca:
            # TODO discuss: does this make sense? cant i just take the depth?
            if getIC(concept, tree, ic_mode) > getIC(lca, tree, ic_mode):
                lca = concept
    return lca

def getAncestors(concept: str, tree: Tree):
    """Return the ancestors of a concept in a given tree"""
    if concept == tree.root:
        return set()
    ancestors = []
    parent: Node = tree.parent(concept)
    while tree.depth(parent.identifier) >= 0:
        ancestors.append(parent.identifier)
        if parent.is_root():
            return set(ancestors)
        parent = tree.parent(parent.identifier)

    return set(ancestors)

def getShortestPath(concept1: str, concept2: str, depth_lca: int, tree: Tree):
    depth_concept1 = tree.level(concept1)
    depth_concept2 = tree.level(concept2)
    return depth_concept1 + depth_concept2 - 2*depth_lca

def getCS(concept1: str, concept2: str, tree: Tree, depth: int,ic_mode: str,cs_mode: str):
    """Returns concept similarity of two concepts based on CS-algorithms from https://doi.org/10.1186/s12911-019-0807-y"""
    
    global cs_table
    
    if concept1 == 'nan' or concept2 == 'nan':
        return 0
    if not isinstance(concept1, str) or not isinstance(concept2,str):
        if np.isnan(concept1) or np.isnan(concept2):
            return 0    
    if cs_table.get((concept1,concept2)):
        return cs_table[(concept1,concept2)]
    if cs_table.get((concept2,concept1)):
        return cs_table[(concept2,concept1)]
    
    if concept1 == concept2:
        if cs_mode == 'wu_palmer' or cs_mode == 'simple_wu_palmer':
            cs = 1.0
            cs_table[(concept1,concept2)] = cs
            cs_table[(concept2,concept1)] = cs
            return cs
        elif cs_mode == 'path_based':
            cs = 0.0
            cs_table[(concept1,concept2)] = cs
            cs_table[(concept2,concept1)] = cs
            return cs
    
    lca = getLCA(concept1, concept2, tree, ic_mode)
    ic_lca = getIC(lca, tree,ic_mode)
    ic_1 = getIC(concept1,tree,ic_mode)
    ic_2 = getIC(concept2,tree,ic_mode)
    
    try:

        if cs_mode == 'path_based':
            cs = cs_algorithms.getPathBasedDist(concept1,concept2,tree,depth)
        elif cs_mode == 'wu_palmer':
            cs = cs_algorithms.getCSWuPalmer(ic_1,ic_2,ic_lca) 
        elif cs_mode == 'li':
            cs = cs_algorithms.getCSLi(ic_1,ic_2,ic_lca)
        elif cs_mode == 'simple_wu_palmer':
            cs = cs_algorithms.getCSSimpleWuPalmer(ic_lca,depth)
        elif cs_mode == 'leacock_chodorow':
            cs = cs_algorithms.getCSLeacockChodorow(ic_1,ic_2,ic_lca,ic_mode,tree,depth)
        elif cs_mode == 'nguyen_almubaid':
            cs = cs_algorithms.getCSNguyenAlMubaid(concept1, concept2, lca, tree, depth) 
        elif cs_mode == 'batet':
            cs = cs_algorithms.getCSBatet(concept1, concept2, tree)        
        else:
         raise ValueError('Unsupported CS-mode: ',cs_mode)

        cs_table[(concept1,concept2)] = cs
        cs_table[(concept2,concept1)] = cs
        return cs
    except ValueError as err:
        print(err.args)
        sys.exit()
    
    ###### ONLY PYTHON >= 3.10.0 #####
    # match cs_mode:
    #     # CS1
    #     case 'binary':
    #         return int(concept1==concept2)
    #     # CS2
    #     case 'wu_palmer':
    #         return 1 - (2*ic_lca)/(ic_1+ic_2) 
    #     # CS 3
    #     case 'li':
    #         return 1 - math.exp( 0.2*(ic_1 + ic_2 - 2*ic_lca) )*
    #     # CS4
    #     case 'simple_wu_palmer':
    #         return (depth - ic_lca)/(depth - 1)

def getSetSim(concepts_1: set, concepts_2: set, setsim_mode: str, tree: Tree, cs_mode: str, ic_mode: str) -> float:
    global depth
    try:
        if len(concepts_1) == 1 and len(concepts_2) == 1:
            return getCS(list(concepts_1)[0],list(concepts_2)[0],tree,depth,ic_mode,cs_mode)
        if len(concepts_1) != 0 and len(concepts_2) != 0:
            if 'nan' in concepts_1 or 'nan' in concepts_2:
                return 0
            if setsim_mode == 'jaccard':
                return setsim_algorithms.getJaccardSetSim(concepts_1, concepts_2)
            elif setsim_mode == 'dice':
                return setsim_algorithms.getDiceSetSim(concepts_1, concepts_2)
            elif setsim_mode == 'cosine':
                return setsim_algorithms.getCosineSetSim(concepts_1, concepts_2)
            elif setsim_mode == 'overlap':
                return setsim_algorithms.getOverlapSetSim(concepts_1, concepts_2)
            elif setsim_mode == 'mean_cs':
                return setsim_algorithms.getMeanCSSetSim(concepts_1, concepts_2, tree, cs_mode, ic_mode)
            elif setsim_mode == 'hierarchical':
                return setsim_algorithms.getHierachicalDistSetSim(concepts_1, concepts_2, tree, cs_mode, ic_mode)
            elif setsim_mode == 'bipartite_matching':
                return setsim_algorithms.getWeightedBipartiteMatchingSim(concepts_1,concepts_2,tree,ic_mode,cs_mode)
            else:
                raise ValueError("Unsupported setsim algorithm: ", setsim_mode)
        else:
            raise ValueError('Empty Concept Set(s)')
    except ValueError as err:
        print(err.args)
        sys.exit()

def getAllConcepts(tree: Tree):
    all_concepts = []
    for node in tree.all_nodes():
        all_concepts.append(node.identifier)
    all_concepts.remove(0)
    return all_concepts

def getDistMatrix(concepts: list, tree: Tree, worker_index, max_workers,ic_mode,cs_mode):
    """
    Function for the parallelized processes. \n 
    Computes the part of the (absolute) distance matrix of the given concepts, 
    that corresponds to the worker index of the calling process.
    """
    depth = tree.depth()
    length = len(concepts)
    start = getStart(worker_index, max_workers, length)
    stop = getStop(worker_index, max_workers, length)
    dist_matrix = np.zeros(shape=(stop-start, length))
    i = 0
    for concept1 in concepts[start:stop]: 
        concept1_index = concepts.index(concept1)
        for concept2 in concepts[concept1_index:]:
            cs = getCS(concept1, concept2, tree,depth,ic_mode,cs_mode)
            # safe CS values in matrix (only upper triangular)
            dist_matrix[i, concepts.index(concept2)] = cs
        i+=1
    return dist_matrix, worker_index

def getStop(worker_index, max_workers, length):
    """Returns logarithmically spaced stop index"""
    if worker_index == max_workers:
        return length
    return getStart(worker_index + 1, max_workers, length)

def getStart(worker_index, max_workers, length):
    '''Returns logarithmically spaced start index'''
    logspace = getSpacing(max_workers, length)
    return math.ceil(logspace[worker_index-1])

def getSpacing(max_workers, length):
    """Returns spacing for the concept list."""
    logspace =  length/10*np.logspace(start=-1,stop=1,num=max_workers, endpoint=True)
    # remove offset
    logspace = logspace - logspace[0] 
    return logspace

def getDistMatrixWrapper(p):
    """Wrapper for the parallel-process-function"""
    return getDistMatrix(*p)

def getMDSMatrix(dist_matrix: ndarray) -> DataFrame:
    """Computes multi-dimensionally-scaled two-dimensional concept-coordinates based on a pairwise-distance-matrix"""
    # use MDS to compute the relative distances of the distinct concepts
    embedding = MDS(n_components=2)
    dist_matrix_transformed = embedding.fit_transform(dist_matrix)

    df_dist_matrix = pd.DataFrame(dist_matrix_transformed)
    return df_dist_matrix

def mirrorMatrix(dist_matrix:ndarray) -> ndarray:
    """mirrors uppertriangular distance matrix along its diagonal"""
    return dist_matrix + dist_matrix.T - np.diag(np.diag(dist_matrix))

def plotDistMatrix(coordinates: DataFrame, datapoints: list):
    fig, ax = plt.subplots()
    coordinates.plot(0, 1, kind='scatter', ax=ax)
    ax.set_ylabel('')
    ax.set_xlabel('')
    for k, v in coordinates.iterrows():
        ax.annotate(datapoints[k], v)

    plt.show(block=True)

def saveConceptDistancesInExcel(df_mds_coordinates: DataFrame, concepts: list):
    """Saves pairwise concept-distances to excel."""
    array = df_mds_coordinates.to_numpy()
    dm = distance_matrix(array,array)
    df = pd.DataFrame(dm)
    df.to_excel('concept_distances.xlsx')

def getRandomConcepts(concept_cnt: int,tree: treelib.Tree) -> set:
    """Returns list with concept_cnt random concepts from the given taxonomy tree."""
    nodes: list[Node]
    nodes = random.sample(tree.all_nodes(),concept_cnt)
    concepts = [x.identifier for x in nodes]
    return set(concepts)

def getConceptCount(tree: treelib.Tree) -> int:
    """Returns the number of concepts in a taxonomy."""
    return len(tree.leaves())    

def setMaxIC(tree: Tree, ic_mode: str):
    max_ic = 0
    for node in tree.all_nodes():
        concept = node.identifier
        ic = getIC(concept,tree,ic_mode)
        if ic > max_ic:
            max_ic = ic
    tree.create_node('max_ic','max_ic', data=max_ic,parent=0)
    return

def getMaxIC(tree: Tree, ic_mode: str) -> float:
    if not tree.contains('max_ic'):
            setMaxIC(tree, ic_mode)
    return tree.get_node('max_ic').data

def setMaxDistOrSim(tree: Tree, ic_mode: str, cs_mode: str):
    max_dist_or_sim = 0
    depth = tree.depth()
    for node1 in tree.all_nodes():
        concept1 = node1.identifier
        for node2 in tree.all_nodes():
            concept2 = node2.identifier
            cs = getCS(concept1,concept2,tree,depth,ic_mode,cs_mode)
            if cs > max_dist_or_sim:
                max_dist_or_sim = cs
    tree.create_node('max_dist_or_sim','max_dist_or_sim', data=max_dist_or_sim,parent=0)
    return

def getMaxDistOrSim(tree: Tree, ic_mode: str, cs_mode: str) -> float:
    if not tree.contains('max_dist_or_sim'):
            setMaxDistOrSim(tree, ic_mode, cs_mode)
    return tree.get_node('max_dist_or_sim').data

def getCSMatrix(concepts_1: list, concepts_2: list, tree: Tree, ic_mode, cs_mode) -> ndarray:
    """ Returns CS matrix for given concept sets. """
    cs_matrix = np.zeros(shape=(len(concepts_1),len(concepts_2)))
    depth = tree.depth()

    for concept1 in concepts_1:
        c1_index = concepts_1.index(concept1)
        for concept2 in concepts_2:
            c2_index = concepts_2.index(concept2)
            cs_matrix[c1_index,c2_index] = getCS(concept1,concept2,tree,depth,ic_mode,cs_mode)
            
    return cs_matrix

def normalize(matrix: ndarray) -> ndarray:
    if np.max(matrix) == 0:
        return matrix
    return matrix/np.max(matrix)

def getSetDistMatrixWrapper(p):
    """Wrapper for the parallel-process-function"""
    return getSetDistMatrix(*p)

def getSetDistMatrix(sets: list, tree: Tree, worker_index, max_workers,ic_mode,cs_mode,setsim_mode) -> ndarray:
    """
    Function for the parallelized processes. \n 
    Computes the part of the (absolute) distance matrix of the given conceptsets, 
    that corresponds to the worker index of the calling process.
    """
    length = len(sets)
    start = getStart(worker_index, max_workers, length)
    stop = getStop(worker_index, max_workers, length)
    dist_matrix = np.zeros(shape=(stop-start, length))
    i = 0
    for set1 in sets[start:stop]: 
        set1_index = sets.index(set1)
        for set2 in sets[set1_index:]:
            setsim = getSetSim(concepts_1=set1, concepts_2=set2, tree=tree,ic_mode=ic_mode,cs_mode=cs_mode,setsim_mode=setsim_mode)
            # safe SetSim values in matrix (only upper triangular)
            dist_matrix[i, sets.index(set2)] = setsim
        i+=1
    return dist_matrix, worker_index

def getScaledSetSim(setSim, len1, len2):
    ''' Used to scale the set-similarities to account for differences in set-sizes that might impair the accuracy of the calculations. '''
    setDiff = abs(len1 - len2)
    maxSim = min(len1,len2)
    if maxSim == 0:
        return setSim
    return setSim/(maxSim + math.log(setDiff + 1))

def getModifierLabel(root: ET.Element, modifier: str, mod_code: str) -> str:
    for mod_class in root.iter('ModifierClass'):
        if mod_class.get('code') == mod_code and mod_class.get('modifier') == modifier:
            return mod_class.find('Rubric').find('Label').text