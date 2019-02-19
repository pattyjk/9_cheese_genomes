## DESeq 2 analysis to identify 'upregulated' VOCs
```
#R v. 3.4.4

#install DESeq2
#install.packages("htmltools")
#library(htmltools)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

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

#write results to a table
results_de<-results(dds)
results_table<-as.data.frame(results_de)
results_table$Compound<-row.names(results_table)

#write data to data frame
write.table(results_table, 'deseq_results.txt', row.names=F, quote=F, sep='\t')

#identify which are singificant
deseq_sig<-which(results_table$padj<0.05)

#extract VOC of interest
results_table_sig<-results_table[deseq_sig,]

#plot volcano plot
ggplot(results_table, aes(log2FoldChange, pvalue))+
  geom_point(aes(size=2), colour='grey', show.legend = F)+
  geom_point(data=results_table_sig, aes(log2FoldChange, pvalue, size=2), colour='red',show.legend = F)+
  scale_y_log10()+
  theme(axis.text=element_text(size=15, colour="black"), axis.title=element_text(size=20, colour="black"), strip.text = element_text(size = 12))+
  theme_bw()+
  ylab("Log10 p-value")+
  xlab("Log2 fold change")

#write significant VOC to file
write.table(results_table_sig, 'deseq_results_significant.txt', row.names=F, quote=F, sep='\t')
```
