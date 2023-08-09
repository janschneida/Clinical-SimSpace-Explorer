import src.calc_utils as cu
import src.data_utils as du
import pandas as pd
import numpy as np
import csv
#paac_msk_jco_2023, acbc_mskcc_2015,acyc_jhu_2016

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
    # bcc has 16667 genes -> should be enough I guess LOL
    mut_dict = du.getMutationProfileDictFromStudy('bcc_unige_2016')#paac_msk_jco_2023, acbc_mskcc_2015,acyc_jhu_2016
    # safe_uniprot_lookup(mut_dict)
    # safe_oncoKB_lookup()
    # safe_panther_terms(mut_dict)
    # safe_pygosemsim_terms(mut_dict)
    # safe_cbio_study('acbc_mskcc_2015')
    # du.__initOncoKBLookUp__()
    # genes = du.oncoKB_lookup
    genes = du.getUniqueMutations(mut_dict)
    term_dict = du.__initPantherDBtermDict__(genes[:10])
    # NOTE we can right to new line easily withour overriden eversrhing
    # safe_gene_sim_look_up(genes[:1000],term_dict)
    # safe_cbio_study('blca_nmibc_2017')
    # safe_cbio_study('chol_msk_2018')
    # safe_cbio_study('breast_alpelisib_2020')
    # safe_cbio_study('brca_jup_msk_2020')
    # safe_cbio_study('brca_mapk_hp_msk_2021')
    # safe_cbio_study('crc_nigerian_2020')
    # safe_cbio_study('crc_dd_2022')
    # safe_cbio_study()
    # safe_cbio_study()
    