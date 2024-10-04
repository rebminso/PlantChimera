#Rscript

# Load necessary libraries
library(GenomicFeatures)
library(biomaRt)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
gtf_file <- args[1]
fasta_file <- args[2]
output_file <- args[3]

# Read GTF file
gtf_data <- read.csv(gtf_file, sep = '\t')
head(gtf_data)

# Read FASTA file
con <- file(fasta_file, "r")
mydata <- readLines(con)
close(con)

# Extract headers from FASTA file
txStart <- unlist(lapply(mydata, function(x) startsWith(x, ">")))
txStartID <- which(txStart)
txInfo <- mydata[txStartID]

# Parse transcript information
chrInfo <- unlist(lapply(txInfo, function(x) { 
  y <- unlist(strsplit(x, " "))
  z <- unlist(strsplit(y[3], ":"))
  return(z[3])
}))

cdnaType <- unlist(lapply(txInfo, function(x) { 
  y <- unlist(strsplit(x, " "))
  z <- unlist(strsplit(y[grep("cdna", y)], ":"))
  return(z[length(z)])
}))

geneType <- unlist(lapply(txInfo, function(x) { 
  y <- unlist(strsplit(x, " "))
  z <- unlist(strsplit(y[grep("gene_biotype", y)], ":"))
  return(z[length(z)])
}))

txType <- unlist(lapply(txInfo, function(x) { 
  y <- unlist(strsplit(x, " "))
  z <- unlist(strsplit(y[grep("transcript_biotype", y)], ":"))
  return(z[length(z)])
}))

txName <- unlist(lapply(txInfo, function(x) { 
  y <- unlist(strsplit(x, " "))
  return(substring(y[1], 2))
}))

geneName <- unlist(lapply(txInfo, function(x) { 
  y <- unlist(strsplit(x, " "))
  z <- unlist(strsplit(y[grep("^gene:", y)], ":"))
  return(z[length(z)])
}))

txName <- unlist(lapply(txName, function(x) unlist(strsplit(x, "\\."))[1]))
geneName <- unlist(lapply(geneName, function(x) unlist(strsplit(x, "\\."))[1]))

txAnno <- cbind(txName, cdnaType, chrInfo, geneType, txType, geneName)
geneAnno <- txAnno[!duplicated(as.character(txAnno[, 6])), ]

##### Step 2: Using BioMart to get gene paralog and other information
### replace "athaliana_eg_gene" with youe ensembl plant dataset name
ensembl <- useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org")
geneParalog <- getBM(attributes = c('ensembl_gene_id', "athaliana_eg_paralog_ensembl_gene"), 
                     filters = 'ensembl_gene_id', values = geneAnno[, 6], mart = ensembl)

pick <- !is.na(geneParalog[, 2]) & geneParalog[, 2] != ""
geneParalog <- geneParalog[pick, drop = FALSE, ]

# Rename the columns of the geneParalog data frame
colnames(geneParalog) <- c("GeneID", "ParalogGeneID")

write.table(geneParalog, file = output_file, sep = '\t', quote = FALSE, row.names = FALSE)
