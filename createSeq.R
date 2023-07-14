#find which are used for sequencing
all_files <- basename(list.files(path = '~', full.names = TRUE))
value <- "gene_cds"
selected <- all_files[which(grepl(value, all_files))]

library(seqinr)
library(Biostrings)
library(coRdon)
library(dplyr)
library(stringr)


'multiple files'
if (length(selected) > 1) {
  seqs <- c()
  final <- data.frame( sequence = character(0))
  for (items in selected) {
    file_path = paste0("~/",items)
    set = unzip(file_path)[2]
    reading = read.fasta(set)
    seqs <- append(seqs, reading)
  }
  
  gene <- c()
  for ( i in names(seqs)) {
    string <- attr(seqs[[i]], "Annot")
    split <- unlist(strsplit(string, " "))
    gene <- append(gene,split[2])
  }
  
  for ( x in unique(gene)) {
    check <- which(gene == x)
    print(paste("Loading...", x))
    exons2 <- toupper(paste(unlist(seqs[check]), collapse = ""))
    final[x,] = DNAStringSet(exons2)
    check <- c()
  }
  
  finalseq <- DNAStringSet(as.matrix(final))
  codons.df <- as.data.frame(codonTable(finalseq)@counts)
  
  write.csv(final, file = "entrez_seq.csv") 
  
  
  
  
'one file'
} else {
  file_path = paste0("~/", selected)
  set <- unzip(file_path)[2]
  reading = read.fasta(set) #read.fasta(unzip("zipped.zip))
  final <- data.frame( sequence = character(0))
  
  gene <- c()
  for ( i in names(reading)) {
    string <- attr(reading[[i]], "Annot")
    split <- unlist(strsplit(string, " "))
    gene <- append(gene,split[2])
  }
  for ( x in unique(gene)) {
    check <- which(gene == x)
    print(paste("Loading...", x))
    exons2 <- toupper(paste(unlist(reading[check]), collapse = ""))
    final[x,] = DNAStringSet(exons2)
    check <- c()
  }
  
  finalseq <- DNAStringSet(as.matrix(final))
  codons.df <- as.data.frame(codonTable(finalseq)@counts)
  
  write.csv(final, file = "entrez_seq.csv") 
}
