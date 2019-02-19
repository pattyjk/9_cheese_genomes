## DESeq 2 analysis to identify 'upregulated' VOCs
```
#R v. 3.4.4

#install DESeq2
install.packages("htmltools")
library(htmltools)
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

#read in packages
library( "DESeq2" )
library(ggplot2)
#DeSeq2 1.22.2

#read in data and metadata
table<-read.delim("/Users/patty/OneDrive/Desktop/cheese_volitiles.txt", header=T, row.names = 1)
meta<-read.delim("/Users/patty/OneDrive/Desktop/cheese_volitiles_meta.txt", header=T)

#replace 'NA's' with zeroes
table[is.na(table)] <- 0

#change numbers to intergers by rounding



table<-round_df(table, -1)
View(table)

table<-as.integer(table)

#convert to DESeq format
meta_deseq<-DESeqDataSetFromMatrix(countData=table, colData=meta, design = ~Community)

#run deseq pipeline
dds<-DESeq(meta_deseq)
results_de<-results(dds)

#plot data
plotMA(results_de)
```
