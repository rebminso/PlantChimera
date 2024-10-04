#Rscript

# Load necessary libraries
library(GenomicFeatures)
library(biomaRt)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
gtf_file <- args[1]
fasta_file <- args[2]
output_file <- args[3]

# Set working directory (optional, if needed)
# setwd("/home/garima/Downloads/fusion/FUSION_TOOL/outputs/FuSeq/ArabidopsisThaliana")

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
ensembl <- useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org")
geneParalog <- getBM(attributes = c('ensembl_gene_id', "athaliana_eg_paralog_ensembl_gene"), 
                     filters = 'ensembl_gene_id', values = geneAnno[, 6], mart = ensembl)

pick <- !is.na(geneParalog[, 2]) & geneParalog[, 2] != ""
geneParalog <- geneParalog[pick, drop = FALSE, ]

# Rename the columns of the geneParalog data frame
colnames(geneParalog) <- c("GeneID", "ParalogGeneID")

write.table(geneParalog, file = output_file, sep = '\t', quote = FALSE, row.names = FALSE)


# setwd("/home/garima/Downloads/fusion/FUSION_TOOL/outputs/FuSeq/ArabidopsisThaliana")
# getwd()

# gtfFile=read.csv('at.gtf',sep = '\t')
# head(gtfFile)


# #gtfSqliteFn=args[2];
# library(GenomicFeatures);
# # gtfTxdb <- makeTxDbFromGFF(file=gtfFile,
# #                            format="gtf",
# #                            dataSource=paste("Link to the source",sep=""),
# #                            organism="Arabidopsis thaliana")
# # 
# # 
# #help(makeTxDbFromGFF)
# con <- file("Arabidopsis_thaliana.TAIR10.cdna.all.fa", "r")
# con
# head(con)
# mydata= readLines("Arabidopsis_thaliana.TAIR10.cdna.all.fa")
# head(mydata) 
# # extracting headers from fasta file
# txStart=unlist(lapply(mydata, function(x) startsWith(x,">")))
# head(txStart)
# txStartID=which(txStart)
# head(txStartID)
# txInfo=mydata[txStartID]
# head(txInfo)


# chrInfo=unlist(lapply(txInfo, function(x) { 
#   y=unlist(strsplit(x," "))
#   z=unlist(strsplit(y[3],":"))
#   return(z[3])
# }))
# head(chrInfo)

# cdnaType=unlist(lapply(txInfo, function(x) { 
#   y=unlist(strsplit(x," "))
#   z=unlist(strsplit(y[grep("cdna",y)],":"))
#   return(z[length(z)])
# }))
# head(cdnaType)
# geneType=unlist(lapply(txInfo, function(x) { 
#   y=unlist(strsplit(x," "))
#   z=unlist(strsplit(y[grep("gene_biotype",y)],":"))
#   return(z[length(z)])
# }))
# head(geneType)
# txType=unlist(lapply(txInfo, function(x) { 
#   y=unlist(strsplit(x," "))
#   z=unlist(strsplit(y[grep("transcript_biotype",y)],":"))
#   return(z[length(z)])
# }))
# head(txType)
# txName=unlist(lapply(txInfo, function(x) { 
#   y=unlist(strsplit(x," "))
#   return(substring(y[1],2))
# }))
# head(txName)
# geneName=unlist(lapply(txInfo, function(x) { 
#   y=unlist(strsplit(x," "))
#   z=unlist(strsplit(y[grep("^gene:",y)],":"))
#   return(z[length(z)])
# }))
# head(geneName)

# txName=unlist(lapply(txName, function(x) unlist(strsplit(x,"\\."))[1]))
# head(txName)
# geneName=unlist(lapply(geneName, function(x) unlist(strsplit(x,"\\."))[1]))
# head(geneName)

# txAnno=cbind(txName,cdnaType,chrInfo,geneType,txType,geneName)
# head(txAnno)
# geneAnno=txAnno
# head(geneAnno)
# dupID=duplicated(as.character(geneAnno[,6]))
# geneAnno=geneAnno[!dupID,]
# dim(geneAnno)

# geneAnno

# ##### Step 2) using biomart to get gene paralog and other information
# ### get paralog
# # NOTE: 
# # NEED TO CHECK THE ATTRIBUTE in biomart data that containing GENE PARALOG INFORMATION
# # USUALLY it is with the name format as "SPECIES_paralog_ensembl_gene", for example hsapiens_paralog_ensembl_gene, athaliana_eg_paralog_ensembl_gene, etc
# library(biomaRt)
# #geneParalog=matrix("N.A",1,2)
# #colnames(geneParalog)=c('ensembl_gene_id',"SPECIES_paralog_ensembl_gene")
# ##example for Arabidopsis thaliana


# ensembl = useMart(biomart = "plants_mart",dataset="athaliana_eg_gene", host = "plants.ensembl.org")
# geneParalog <- getBM(attributes=c('ensembl_gene_id',"athaliana_eg_paralog_ensembl_gene"),filters = 'ensembl_gene_id', values = geneAnno[,6], mart = ensembl)
# head(geneParalog)

# pick=!is.na(geneParalog[,2])
# pick
# geneParalog=geneParalog[pick,drop=FALSE,]
# pick=geneParalog[,2]!=""
# geneParalog=geneParalog[pick,drop=FALSE,]
# head(geneParalog)


# str(geneParalog)

# write.table(geneParalog,file = '/home/garima/Documents/Fusion/redo_fusion_June/merged_arab/GeneParalogue_jul11.txt', sep = '\t', quote = FALSE, row.names = FALSE)

