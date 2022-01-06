#import libraries 
library(devtools)
library(dplyr)
library(affy)
library(Biobase)
library(GEOquery)
library(arrayQualityMetrics)
library(BiocManager)
library(oligo)
library(ggplot2)
library(data.table)

#load files
setwd("~/Downloads/TAC_Files/CEL")
CelFiles = list.celfiles("~/Downloads/TAC_Files/CEL", 
                   full.names = TRUE)

#read raw CEL files
rawData <- read.celfiles(CelFiles)

#apply RMA normalization function 
esetAB <- rma(rawData)

#write normalized expression data to txt to add id var in .txt file (manually)
write.exprs(esetAB,"rma_cel.txt")

#load gene ids & expression data
gene_ids <- read.delim("~/Downloads/TAC_Files/GPL23038_bioProcess.an.txt", check.names = FALSE)
data <- read.delim("~/Downloads/TAC_Files/CEL/rma_cel.txt", check.names = FALSE)

#merge gene ids and expression data
all <- left_join(gene_ids, data, by ="ID")

#uppercase genes
all$GeneSymbols <- toupper(all$GeneSymbols)

#remove genes without GeneSymbols and psuedo genes
all <- all %>%
    filter(GeneSymbols != 'NA') %>%
    filter(GeneNames == 'psuedo')  #%>%
    select(-ID,-GOTerms,-GemmaIDs,-NCBIids)

#write final dataset to .csv file
write.csv(all,"annotated_cel.csv")
