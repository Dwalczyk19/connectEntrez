#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 12:19:06 2023

@author: dwalz
#handle = Entrez.esearch(db = "gene", term = ("(Homo sapiens[Organism] OR Homo sapiens[All Fields]) AND CD248[GENE]"))
#Entrez.read(handle)
#(("Homo sapiens"[Organism] OR Homo sapiens[All Fields]) AND CD248[All Fields]) AND alive[prop]
"""

#remeber python starts at 0 not 1
import pandas as pd 
import Bio
from Bio import Entrez
from Bio import SeqIO

read = pd.read_excel("/home/dwalz/R/data/Human_genes_with_entrez_IDs060523.xlsx")
read = pd.DataFrame(read)
entrez_id = read["NEW-Entrez-ID"]
gene_id = read["NEW-Gene-ID"]

#for google 
from google.colab import drive
drive.mount("/content/drive")

import os 
path = "/content/drive/MyDrive/rna institue data/data/Human_genes_with_entrez_IDs060523.xlsx"
df = pd.read_excel(path)


Entrez.email = "walczd3@rpi.edu"

def get(entrez_id): 
  handle = Entrez.efetch(db = "gene", id = "57124", rettype = "fasta", retmode = "text")
  coin = handle.read().strip()
  area = coin.split("\n")[4:6]
  location = area[0].split("; ")[1].split(": ")[1] #first number references that its on chromosome 'x'
  #get location and check if chromosome has already been mapped to save data
  #apply map(func, list) to apply function to data 
  #once sequence is found add it to dataframe with gene in col 1 & sequence in col 2 (or do this separately and just get sequence for now)





