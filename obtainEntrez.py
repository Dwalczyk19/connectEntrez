#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 12:19:06 2023

@author: dwalz
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

test = entrez_id[0]
Entrez.email = "walczd3@rpi.edu"
handle = Entrez.esearch(db = "gene", term = ("(Homo sapiens[Organism] OR Homo sapiens[All Fields]) AND CD248[GENE]"))
Entrez.read(handle)
#(("Homo sapiens"[Organism] OR Homo sapiens[All Fields]) AND CD248[All Fields]) AND alive[prop]
