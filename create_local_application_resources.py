import src.calc_utils as cu
import src.data_utils as du
import pandas as pd
import numpy as np
import csv

def safe_panther_terms(mut_dict):
    
    du.__initPantherDBtermDict__(du.getUniqueMutations(mut_dict))
    with open('panther_term_dict.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header row
        header = ['Gene', 'GO Terms']
        writer.writerow(header)

        # Write the gene-GO term associations
        for gene, terms in du.panther_terms.items():
            row = [gene, ' '.join(terms)]
            writer.writerow(row)   
    
def safe_pygosemsim_terms(mut_dict):
    du.__initPyGOSemSimTermDict__(du.getUniqueMutations(mut_dict))
    with open('pygosemsim_term_dict.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header row
        header = ['Gene', 'GO Terms']
        writer.writerow(header)

        # Write the gene-GO term associations
        for gene, terms in du.pygosemsim_terms.items():
            row = [gene, ' '.join(terms)]
            writer.writerow(row)   
            
def safe_uniprot_lookup(mut_dict):
    du.__initUniprotLookUp__(du.getUniqueMutations(mut_dict))
    with open('uniprot_lookup.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header row
        header = ['Uniprot', 'HUGO Gene Name']
        writer.writerow(header)

        for uniprot, hugo_name in du.uniprot_lookup.items():
            row = [uniprot, hugo_name]
            writer.writerow(row)   

def safe_oncoKB_lookup():
    du._initOncoKBLookUp_()
    with open('oncoKB_lookup.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header row
        header = ['Oncogene/TSG Gene']
        writer.writerow(header)

        for gene in du.oncoKB_lookup:
            row = [gene]
            writer.writerow(row)   

def safe_cbio_study(study):
    df = du._getDataFrameFromCBioPortalStudy_(study)
    df.to_csv(study)
    
def safe_gene_sim_look_up(genes,term_dict):
    du._initPyGoGO_()
    cu.__initGeneSimLookUp__(genes,term_dict)
    with open('gene_sim_lookup.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header row
        header = ['Gene Pair','Gene Sim']
        writer.writerow(header)

        for genePair,geneSim in du.gene_sim_look_up.items():
            row = [' '.join(genePair),geneSim]
            writer.writerow(row)


if __name__ == '__main__':
    return
    
