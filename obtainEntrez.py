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


read = pd.read_excel("~/data/Human_genes_with_entrez_IDs060523.xlsx")
read = pd.DataFrame(read)
entrez_id = read["NEW-Entrez-ID"]
gene_id = read["NEW-Gene-ID"]


configuration = ncbi.datasets.openapi.Configuration(
    host = "https://api.ncbi.nlm.nih.gov/datasets/v1"
)

configuration.api_key['ApiKeyAuthHeader'] = 'YOUR API KEY'

Entrez.email = simpledialog.askstring(title = "Entrez Email", prompt = "Please enter your email for NCBI authentication")
    

def get_gene(gene_ids, mcheck):
    
    if mcheck == False :
        n = len(gene_ids)
        L = []
        
        if n > 1000: 
            if n % 1000 == 0: 
                steps =int( n / 1000)
                for i in range(1,steps+1):
                    ct = i * 1000
                    L.append((ct-1000,ct))
                
            else: 
                steps = int(n / 1000)
                rm = n % 1000
                
                if steps == 1: 
                    L.append((0,1000))
                    L.append((1000,1000+rm))
                else: 
                    for i in range(1,steps+1): 
                        ct = i * 1000 
                        L.append((ct-1000,ct))
                        if i == steps: 
                            L.append((ct, ct + rm))
            
            for items in L: 
                num = f'{items[0], items[1]}'
                zipfile_name = "gene_cds" + num + ".zip"
                with DatasetsApiClient() as api_client:
                    
                    gene_api = DatasetsGeneApi(api_client)
                    try:
                        print("Starting download of...{}".format(zipfile_name))
                        gene_dataset_download = gene_api.download_gene_package(
                            gene_ids[int(items[0]):int(items[1])],
                            include_annotation_type=["FASTA_CDS"],
                            _preload_content=False,
                        )
                        
                        with open(zipfile_name, "wb") as f:
                            f.write(gene_dataset_download.data)
                        print("Finished Downloading...")
                    except DatasetsApiException as e:
                        sys.exit(f"Exception when calling GeneApi: {e}\n")
        else: 
            num = "(0-" + str(n) + ")"
            zipfile_name = "gene_cds" + num + ".zip"
            with DatasetsApiClient() as api_client:
                
                gene_api = DatasetsGeneApi(api_client)
                try:
                    print("Starting download of...{}".format(zipfile_name))
                    gene_dataset_download = gene_api.download_gene_package(
                        gene_ids,
                        include_annotation_type=["FASTA_CDS"],
                        _preload_content=False,
                    )
                    
                    with open(zipfile_name, "wb") as f:
                        f.write(gene_dataset_download.data)
                    print("Finished Downloading...")
                except DatasetsApiException as e:
                    sys.exit(f"Exception when calling GeneApi: {e}\n")
                    
    else: 
        zipfile_name = "gene_cds.zip"
        with DatasetsApiClient() as api_client:
            
            gene_api = DatasetsGeneApi(api_client)
            try:
                print("Starting download of...{}".format(zipfile_name))
                gene_dataset_download = gene_api.download_gene_package(
                    gene_ids,
                    include_annotation_type=["FASTA_CDS"],
                    _preload_content=False,
                )
                
                with open(zipfile_name, "wb") as f:
                    f.write(gene_dataset_download.data)
                print("Finished Downloading...")
            except DatasetsApiException as e:
                sys.exit(f"Exception when calling GeneApi: {e}\n")
        
    
get_gene(list(entrez_id[:<ENTER YOUR VALUE>]), mcheck = False)




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
    
#Enter your manual values 
gene_ids = []  
get_gene(gene_ids, mcheck = True)




