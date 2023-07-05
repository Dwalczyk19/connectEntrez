read = pd.read_excel("/home/dwalz/R/data/Human_genes_with_entrez_IDs060523.xlsx")
read = pd.DataFrame(read)
entrez_id = read["NEW-Entrez-ID"]
gene_id = read["NEW-Gene-ID"]

from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GeneApi as DatasetsGeneApi

# Provide your own gene ids as a list of integers

def example_usage_of_api(gene_ids):
    if len(gene_ids) == 0:
        print("Please provide at least one gene-id")
        return

    with DatasetsApiClient() as api_client:
        gene_api = DatasetsGeneApi(api_client)

        # Get just metadata
        try:
            gene_reply = gene_api.gene_metadata_by_id(gene_ids)
            for gene in gene_reply.genes:
                print(gene.gene.gene_id)
        except DatasetsApiException as e:
            print(f"Exception when calling GeneApi: {e}\n")

        # Or, download a data package with FASTA files
        try:
            print("Begin download of data package ...")
            gene_ds_download = gene_api.download_gene_package(
                gene_ids, include_annotation_type=["FASTA_GENE"], _preload_content=False
            )
            gene_reply = gene_api.gene_metadata_by_id(gene_ids)
            zipfile_name = "sequence_info.zip"

            with open(zipfile_name, "wb") as f:
                f.write(gene_ds_download.data)
            print(f"Download completed -- see {zipfile_name}")

        except DatasetsApiException as e:
            print(f"Exception when calling GeneApi: {e}\n")
            

            
example_usage_of_api(list(entrez_id))




