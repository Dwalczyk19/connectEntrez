
import ncbi.datasets.openapi    
from ncbi.datasets.openapi.api import gene_api
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GeneApi as DatasetsGeneApi
from ncbi.datasets.package import dataset
from Bio import Entrez, SeqIO
from numpy import array as arr
import pandas as pd
import regex as re 


# Provide your own gene ids as a list of integers
read = pd.read_excel("/home/dwalz/R/data/Human_genes_with_entrez_IDs060523.xlsx")
read = pd.DataFrame(read)
entrez_id = read["NEW-Entrez-ID"]
gene_id = read["NEW-Gene-ID"]


Entrez.email = "walczd3@rpi.edu"
configuration = ncbi.datasets.openapi.Configuration(
    host = "https://api.ncbi.nlm.nih.gov/datasets/v1"
)
configuration.api_key['ApiKeyAuthHeader'] = "9cbb475748fcce2a676126874b4fc0616f08"
zipfile_name = "sequence_info.zip"

def example_usage_of_api(gene_ids):
    #if len(gene_ids) == 0:
     #   print("Please provide at least one gene-id")
      #  return

    with DatasetsApiClient() as api_client:
        gene_api = DatasetsGeneApi(api_client)

        try:
            print("Begin download of data package ...")
            gene_ds_download = gene_api.download_gene_package(
                gene_ids, include_annotation_type=["FASTA_GENE"], _preload_content=False   #FASTA_PROTEIN, FASTA_RNA, etc.
            )
            data = gene_ds_download.data
            #store = "ncbi_dataset/data/data_report.jsonl".join(str(data).split("ncbi_dataset/data/data_report.jsonl")[:2])
            #print(store.decode("utf-8"))
            
            with open(zipfile_name, "wb") as f:
                f.write(data)
            print(f"Download completed -- see {zipfile_name}")
            

        except DatasetsApiException as e:
            print(f"Exception when calling GeneApi: {e}\n")
            

n = int(input("How many genes would you like to access, press 0 to quit --> "))
if n > 0: 
    example_usage_of_api(list(entrez_id)[:3])  #specify index before running, *1050 but use 1000 for simplicity, 12 times to retrieve the entire file. 11910 rows (nearly 12k)
elif n == 0: 
    print("Finished...")
else: 
    n = int(input("How many genes would you like to access, press 0 to quit --> "))
    


