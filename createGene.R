library(seqinr)
library(Biostrings)
library(coRdon)
library(dplyr)
library(stringr)

file_path = "~/sequence_info.zip"
set <- unzip(file_path)[2]
reading = read.fasta(set) #read.fasta(unzip("zipped.zip))

#goal is to merge exons together from each 
reading[c(1,2)]
gene <- c()
for ( i in names(reading)) {
  string <- attr(reading[[i]], "Annot")
  split <- unlist(strsplit(string, " "))
  gene <- append(gene,split[2])
}

seq <- ""

#seq <- gsub("[\r\n]", "", seq) #removes \n
#which(unlist(strsplit(seq, split = "")) != unlist(strsplit(unL, split = "")))
unL <- toupper(paste(unlist(reading[26]), collapse = ""))
final <- data.frame( sequence = character(0))

for ( x in unique(gene)) {
  check <- which(gene == x)
  print(paste("Loading...", x))
  exons2 <- toupper(paste(unlist(reading[check[1]]), collapse = ""))
  final[x,] = exons2
  check <- c()
}


write.csv(final, file = "entrez_seq.csv")


