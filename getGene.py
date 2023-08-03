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
import jsonlines
from zipfile import ZipFile
import os



            
def download(zipfile_name, items, gene_ids, m):
    
    #m = 1 = > 1000;  2 = < 1000;  3 = manual
    
    if m == 1: 
    
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
    elif m == 2 or m == 3: 
        with DatasetsApiClient() as api_client:
            
            gene_api = DatasetsGeneApi(api_client)
            try:
                print("Starting download of...{}".format(zipfile_name))
                gene_dataset_download = gene_api.download_gene_package(
                    gene_ids,
                    include_annotation_type=["FASTA_CDS", "FASTA_GENE", "FASTA_RNA"],
                    _preload_content=False,
                )
                
                with open(zipfile_name, "wb") as f:
                    f.write(gene_dataset_download.data)
                print("Finished Downloading...")
            except DatasetsApiException as e:
                sys.exit(f"Exception when calling GeneApi: {e}\n")
    elif m == 4: 
        with DatasetsApiClient() as api_client:
            
            
            gene_api = DatasetsGeneApi(api_client)
            gene_dataset_download = gene_api.download_gene_package(
                gene_ids,
                include_annotation_type=["FASTA_GENE", "FASTA_RNA"],
                _preload_content=False,
            )
            
            with open(zipfile_name, "wb") as f:
                f.write(gene_dataset_download.data)
            

    
            
            
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
                download(zipfile_name, items, gene_ids, 1)
                        
        else: 
            num = "(0-" + str(n) + ")"
            zipfile_name = "gene_cds" + num + ".zip"
            download(zipfile_name, None, gene_ids, 2)
                    
    else: 
        zipfile_name = "gene_cds.zip"
        download(zipfile_name, None, gene_ids,3)



def checkAssembly(file, report_file):

    gene_file = {} 

    #file
    with jsonlines.open(report_file) as reader:

    
        for obj in reader.iter(type=dict):     

            try:
                gene_file[obj["geneId"]] = [ obj["transcripts"][0]["cds"]["accessionVersion"]] #length one 
                
            except KeyError as e: 
                e = str(e)
                
                if e == "'cds'":
                    if obj["geneId"] == "55199": #change to iterate over all types of genes (this is psuedogene & id specific)
                        gene_file[obj["geneId"]] = [obj["transcripts"][0]["type"], "FASTA_GENE"]
                        
                    else: 
                        gene_file[obj["geneId"]] = [obj["type"], "CDS"]
                    
                
                elif e == "'transcripts'":
                    gene_file[obj["geneId"]] = [obj["type"], "FASTA_GENE"] #length two with FASTA_GENE at index 1
            
    return gene_file

def geneCheck(gene_fna, rna_fna, home, i):
    
    rnaCheck = set()
    rna_df = pd.DataFrame()
    gene_df = pd.DataFrame()

    if rna_fna != "blank":
        with open(home + "/assembly" + str(i) + "/ncbi_dataset/data/rna.fna", "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                transcript = str(record.id).split(":")[0] #NM/XM, NR/XR
                geneID = str(record.description).split("] [")[1][7:-1] #Gene ID
                geneName = str(record.description).split()[1].strip() #Gene name 
                seq = str(record.seq) #sequence
                length = len(seq) #length
                rnaCheck.add(geneID)
                rna_df = pd.concat([rna_df, pd.DataFrame({"Transcript": [transcript], "ID": geneID, "Name":geneName, "Sequence": seq, "Length": length})], ignore_index=True)
                
        
        rna_df = rna_df.groupby("ID", as_index=False).first()
        

    with open(home + "/assembly" + str(i) + "/ncbi_dataset/data/gene.fna", "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            transcript = str(record.id).split(":")[0] #NM/XM, NR/XR
            geneID = str(record.description).split("] [")[1][7:] #Gene ID
            geneName = str(record.description).split()[1].strip() #Gene name 
            seq = str(record.seq) #sequence
            length = len(seq) #length

            if geneID not in rnaCheck: 
                gene_df = pd.concat([gene_df, pd.DataFrame({"Transcript": [transcript], "ID": geneID, "Name":geneName, "Sequence": seq, "Length": length})], ignore_index=True)
                
            else: 
                continue 

    gene_df = gene_df.groupby("ID", as_index = False).first()
    
    return pd.concat([gene_df, rna_df])
            
            
            
            


def findAssembly(L, home, gene_ids):  #convert files in L to something that can be altered and changed into a sequence. Maybe make an initial step to check the validity of each gene via some datasets function
    #for .jsonl
    
    complete = pd.DataFrame()
    from operator import itemgetter
    if len(L) > 1: 
        final = pd.DataFrame()
        for i in range( len( L )): 
            print(L[i], "============================================")
            with ZipFile(L[i], "r") as zip: 
                direct = zip.namelist()
                report = direct[2]
                cds_file = direct[1]
                zip.extract(report, path = "assembly" + str(i))
                zip.extract(cds_file, path = "assembly" + str(i))
            
            
            new_file = home + "/assembly" + str(i)+ "/ncbi_dataset/data/data_report.jsonl"
            assembly = home + "/assembly" + str(i)+ "/ncbi_dataset/data/"
            
            gene_file = checkAssembly(L[i], new_file)
            
            genesort = list(gene_file.keys())
            genesort.sort()

            gene_dict = {int(i): gene_file[i] for i in genesort} 
            #s = set([(y,x[0]) for y, x in gene_dict.items() if x[0] in ["PSEUDO", "ncRNA"]])
            #blank.append(s)
            
            #get CDS gene first this is a majority of the file. get longest transcript associated with each gene id 
            cds = [y[0] for x,y in gene_dict.items() if len(y) == 1]
            assembly = [x for x,y in gene_dict.items() if (len(y) > 1 and y[1] == "CDS")]
            gene_fna = [x for x,y in gene_dict.items() if (len(y) > 1 and y[1] == "FASTA_GENE")]
            cds_df = pd.DataFrame()
            test_df = pd.DataFrame()
            gene_df = pd.DataFrame()
            cds_path = home + "/assembly" + str(i) + "/ncbi_dataset/data/" + "cds.fna"

            with open(cds_path, "r") as handle: #main part
                final = []
                check = 0
                count = 0
                find_length = []
                for record in SeqIO.parse(handle, "fasta"):
                    
                    transcript = str(record.id).split(":")[0] #NM/XM, NR/XR
                    geneID = str(record.description).split("] [")[1][7:] #Gene ID
                    geneName = str(record.description).split()[1].strip() #Gene name 
                    seq = str(record.seq) #sequence
                    length = len(seq) #length
                    
                    if transcript in cds: 
                        cds_df = pd.concat([cds_df, pd.DataFrame({"Transcript": [transcript], "ID": geneID, "Name":geneName, "Sequence": seq, "Length": length})], ignore_index=True)
                    
                    
                    elif int(geneID) in assembly:
                        check = int(geneID)
                        if len(find_length) == 0:
                            find_length.append(check)
                            find_length.append( (transcript, length, geneName, seq))
                        elif check in find_length:
                            find_length.append((transcript, length, geneName, seq))
                        elif check not in find_length:
                            name = find_length[0]
                            find_length.pop(0)
                            top = max(find_length, key = itemgetter(1))
                            test_df = pd.concat([test_df, pd.DataFrame({"Transcript": [top[0]], "ID":name, "Name": top[2], "Sequence": top[3], "Length":top[1]})], ignore_index=True)
                            find_length = []
                    
                    elif int(geneID) in gene_fna:
                        print(int(geneID))
                        
                merge = pd.concat([test_df, cds_df])
            
            #genes/rna
            if len(gene_fna) > 0: 
                gene_path = home + "/assembly" + str(i) + "/ncbi_dataset/data/" + "GENE.zip" 
                download(gene_path, None, gene_fna, 4)
                #extract gene file (from in this case assembly 1)
                #get the gene file, se SeqIO.read(file, 'fasta') to get fasta sequence and store it in gene dataframe 
                
                with ZipFile(gene_path, "r") as zip:
                    direct = zip.namelist()
                    gene_fna = direct[1]
                    rnaCheck = False
                    if "rna" in direct[2]:
                        rnaCheck = True
                        rna_fna = direct[2]
                        zip.extract(rna_fna, path = home + "/assembly" + str(i))
                    zip.extract(gene_fna, path = home + "/assembly" + str(i))
                
                os.remove(gene_path)
                if rnaCheck == True:
                    gene_rna_final = geneCheck(gene_fna, rna_fna, home, i)
                else: 
                    gene_rna_final = geneCheck(gene_fna, "blank", home, i)
                    
                
            pre_final = pd.concat([merge, gene_rna_final])
            
            complete = pd.concat([complete, pre_final])
            
       
        '''#ncRNA introns for id's like 220158 require rna format or gene format to download: no cds'''
        '''finding a way to revamp the whole process of finding the type and if it lies between 
        certain categories than you can download that specific annotation type would be ideal'''
        
        '''also bring it to a specific function so you don't have to worry about typing it over and over'''
            
    
    else: 
        #this extracts the json file into a new folder called assembly
        with ZipFile(L[0], "r") as zip: 
            direct = zip.namelist()
            report = direct[2]
            cds_file = direct[1]
            zip.extract(report, path = "assembly")
            zip.extract(cds_file, path = "assembly")
            
        new_file = home + "/assembly/ncbi_dataset/data/data_report.jsonl"
        with jsonlines.open(new_file) as reader: 
            narray = [obj["transcripts"][0]["cds"]["accessionVersion"] for obj in reader.iter(type = dict)]

                
        t_df = pd.DataFrame(narray, columns = ["Transcript"])
        t_df.to_csv(home + "/assembly/ncbi_dataset/data/transcript.csv", index = False, header = True)
        
        #getting GRCh38.p14 assembly 
    complete[["ID", "Name", "Sequence"]].to_csv("Human-Entrez-IDs.csv", index = False)
    

def merge(): 
    #delete assembly files after merger 
    pass


                
if __name__ == "__main__": 
    pass
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    