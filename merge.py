import sys
import ncbi.datasets.openapi    
from ncbi.datasets.openapi.api import gene_api
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GeneApi as DatasetsGeneApi
from Bio import Entrez, SeqIO
from numpy import array as arr
import pandas as pd 
from tkinter import simpledialog
import os
import shutil
import random

def find_entrez(gene_name):
    pass #use for finding all unfound/incorrect types of geneids 
    

if __name__ == "__main__":
    home = os.path.expanduser("~")
    os.chdir(home)
    sys.path.append(home + "/connectEntrez") #just to get files
    os.chdir(os.getcwd() + "/connectEntrez")
    from getGene import get_gene
    from getGene import findAssembly


    '''if in conda env I add this'''
    
    read = pd.read_excel("~/data/Human_genes_with_entrez_IDs060523.xlsx")
    read = pd.DataFrame(read)
    entrez_id = read["NEW-Entrez-ID"]
    gene_id = read["NEW-Gene-ID"]
    
    configuration = ncbi.datasets.openapi.Configuration(
        host = "https://api.ncbi.nlm.nih.gov/datasets/v1"
    )
    
    configuration.api_key['ApiKeyAuthHeader'] = '<your api key>'
    Entrez.email = simpledialog.askstring(title = "Entrez Email", prompt = "Please enter your email for NCBI authentication")
    
    '''pick n'''
    L = list(entrez_id[random.sample(range(0,11900),11000)]) #test
    get_gene(pd.unique(entrez_id).tolist(), mcheck = False)
    

    '''for CDS'''
    #unzip files (1)
    home = os.path.expanduser("~/connectEntrez")
    complete = os.listdir(home)
    idx = [complete[complete.index(items)]  for items in complete if "gene_cds" in items]
    finding = findAssembly(idx, home, pd.unique(entrez_id).tolist())
    
    
    
    #psuedo, protein_coding, ncRNA
    #for all PSEUDO & ncRNA get the longest transcript of files in rna.fna
 



    '''manual entry'''         
    gene_ids = []  #put ids in here
    check = 0
    while (check != -1): 
        check = simpledialog.askinteger(title = "Initial Entrez ID", prompt = "\t\t\tEnter Entrez ID(s) for download. Press -1 to quit.\nIf your downloading a large number of IDs place them in a list and then use funciton get_gene")
        gene_ids.append(check)
        if check == -1:
            if len(gene_ids) > 0 and gene_ids[0] != -1: 
                gene_ids.remove(-1) 
                get_gene(gene_ids, mcheck = True)
            else:
                gene_ids = []
                break 
    

    #shutil.rmtree(home + "/assembly") #remvoes directory (/ies), be careful!!!! 
    sys.path.remove(home)  

