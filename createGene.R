library(seqinr)
library(Biostrings)
library(coRdon)
library(dplyr)
library(stringr)

file_path = "/home/dwalz/R/data/ncbi_dataset/data/gene.fna"

reading = read.fasta(file_path)

#goal is to merge exons together from each 
reading[1]

gene <- c()
for ( i in names(reading)) {
  string <- attr(reading[[i]], "Annot")
  split <- unlist(strsplit(string, " "))
  gene <- append(gene,split[2])
}

for ( x in unique(gene)) {
  check <- which(gene == x)
  print(x)
  print(check)
  check <- c()
}




codonTable(DNAStringSet(paste(as.character(reading$`NC_060935.1:c66310697-66308147`), collapse = "")))




# Quick examples
# Using paste() function
v <- c('A','B','C','D')
print(v)
print(paste(v,collapse=''))
print(paste(v))



