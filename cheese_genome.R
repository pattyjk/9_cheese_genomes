##9 cheese pangenomes
setwd("/Users/patty/Dropbox/R/cheese_genomes/")

#read in data frame
gene_pa<-read.delim("gene_presence_absence.Rtab", header=T)
dim(gene_pa)
#18434 by 10

#fix community names
names(gene_pa)<-c("gene", "com_341", "com_738", "com_862", "com_876", "com_900", "com_908", "com_947", "com_962", "com_JB" )


#remove rows with group, aka hypothetical genes
gene_pa<-gene_pa[!grepl("group_", gene_pa$gene),]
dim(gene_pa)
#6724 by 10

#write gene names to a file
write.table(gene_pa$gene, "all_genes.txt", row.names=F, quote=F)

#remove gene column
row.names(gene_pa)<-gene_pa$gene
gene_pa<-gene_pa[,-1]
#gene_pa$gene<-row.names(gene_pa)

#make heatmap no dendrogram
#heatmap(as.matrix(gene_pa), Colv = NA, Rowv = NA, col = c("grey", "red"), scale="none", ylab="Gene")

#sample dendrogram
#heatmap(as.matrix(gene_pa), Rowv = NA, col = c("grey", "red"), scale="none", ylab="Gene")

#identify genes conserved accross all 9 comunities
length(which(rowSums(gene_pa) == 9))
#4451

#write core genes to frame
core_genes<-as.data.frame(which(rowSums(gene_pa) == 9))

#remove all other genes from table
core_table<-subset(gene_pa, row.names(gene_pa) %in% row.names(core_genes))
core_table$genes<-row.names(core_table)

#fix gene names
core_table$genes<-gsub('_[[:digit:]]+', '', core_table$genes)
gene_pa$gene<-gsub('_[[:digit:]]+', '', gene_pa$gene)

#write gene names to a frame and file
core_genes<-core_table$genes
write.table(core_genes, "core_genes.txt", quote=F, row.names = F)

#identify rows that only have a sum of 1 (aka only in one community)
length(which(rowSums(gene_pa) == 1))
#724

unique_gene_list<-as.data.frame(which(rowSums(gene_pa) == 1))
dim(unique_gene_list)
#724 by 1 (721 total unique genes)

#make a table of only unique genes
unique_genes<-subset(gene_pa, row.names(gene_pa) %in% row.names(unique_gene_list))
unique_genes$gene<-row.names(unique_genes)

#remove extra stuff from gene names
unique_genes$gene<-gsub("_\\d+", "", unique_genes$gene)

length(unique(unique_genes$gene))
#615 unique genes 

#remove duplicate gene names
unique_genes<-unique_genes[!duplicated(unique_genes$gene),]
nrow(unique_genes)
#615 unique genes

#write unique genes to file
write.table(unique_genes$gene, "unique_gene_list.txt", quote=F, row.names=F)

#identify genes unique to each community
library(dplyr)
com_341<- unique_genes[unique_genes$com_341 == 1,] %>% select("com_341", "gene")
nrow(com_341)
#70 unique genes
com_738<-unique_genes[unique_genes$com_738 == 1,] %>% select("com_738", "gene")
nrow(com_738)
#87 unique genes
com_862<-unique_genes[unique_genes$com_862 == 1,] %>% select("com_862", "gene")
nrow(com_862)
#38 unique genes
com_876<-unique_genes[unique_genes$com_876 == 1,] %>% select("com_876", "gene")
nrow(com_876)
#52 unique genes
com_900<-unique_genes[unique_genes$com_900 == 1,] %>% select("com_900", "gene")
nrow(com_900)
#85 unique genes
com_908<-unique_genes[unique_genes$com_908 == 1,] %>% select("com_908", "gene")
nrow(com_908)
#62 unique genes
com_947<-unique_genes[unique_genes$com_947 == 1,] %>% select("com_947", "gene")
nrow(com_947)
#66 unique genes
com_962<-unique_genes[unique_genes$com_962 == 1,] %>% select("com_962", "gene")
nrow(com_962)
#42 unique genes
com_jb<-unique_genes[unique_genes$com_JB == 1,]%>% select("com_JB", "gene")
nrow(com_jb)
#113 unique genes


######Read in cheese gene information
setwd("./annotations/")
library(plyr)
cheese<-dir(pattern = "\\.tsv$")

#add file names to list of files
names(cheese) <- basename(cheese)

#read in the data and catenate into a single data frame (file name is .id)
cheese_gen <- ldply(cheese, read.delim)

#rename headers
names(cheese_gen)<-c("genome", "locus_tag", "ftype", "bp", "gene", "EC", "COG", "product")

#write to a files
write.csv(cheese_gen, "cheese_genes.csv", quote=F, row.names=F)

#number of genomes in analysis
length(unique(cheese_gen$genome))
#27

#remove '.tsv' from genomes
cheese_gen$genome<-gsub(".tsv", "", cheese_gen$genome)
setwd("/Users/patty/Dropbox/R/cheese_genomes")

#remove extra characters from genes
cheese_gen$gene<-gsub("_[[:digit:]]+", "", cheese_gen$gene)

#remove rows with no genes
cheese_gen<- cheese_gen[-which(cheese_gen$gene == ""),]

#number of unique genes
length(unique(cheese_gen$gene))
##3090

##Add gene info to community lists
head(core_table)

core_table<-merge(core_table, cheese_gen, by.x="genes", by.y="gene", all.x=T)
dim(core_table)
#11902 by 24

#remove duplicate gene names and write to file
core_table<-core_table[!duplicated(core_table$genes),]
dim(core_table)
#4451 by 17
write.table(core_table, "core_table.txt", row.names=F, quote=F, sep="\t")

#add genes to community unique genes
com_341<-merge(com_341, cheese_gen,  by="gene", all.x=T)
com_738<-merge(com_738, cheese_gen,  by="gene", all.x=T)
com_862<-merge(com_862, cheese_gen,  by="gene", all.x=T)
com_876<-merge(com_876, cheese_gen,  by="gene", all.x=T)
com_900<-merge(com_900, cheese_gen,  by="gene", all.x=T)
com_908<-merge(com_908, cheese_gen,  by="gene", all.x=T)
com_947<-merge(com_947, cheese_gen,  by="gene", all.x=T)
com_962<-merge(com_962, cheese_gen, by="gene", all.x=T)
com_jb<-merge(com_jb, cheese_gen, by="gene", all.x=T)

#remove duplicate gene names
com_341<-com_341[!duplicated(com_341$gene),]
com_738<-com_738[!duplicated(com_738$gene),]
com_862<-com_862[!duplicated(com_862$gene),]
com_876<-com_876[!duplicated(com_876$gene),]
com_900<-com_900[!duplicated(com_900$gene),]
com_908<-com_908[!duplicated(com_908$gene),]
com_947<-com_947[!duplicated(com_947$gene),]
com_962<-com_962[!duplicated(com_962$gene),]
com_jb<-com_jb[!duplicated(com_jb$gene),]


#map genes back to KEGG
library(readr)
full_kegg <- read_delim("full_kegg.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
dim(full_kegg)

com_jb<-merge(com_jb, full_kegg, by.x='gene', by.y='Gene', all.x=T, no.dups=T)
com_341<-merge(com_341, full_kegg, by.x='gene', by.y='Gene', all.x=T, no.dups=T)
com_738<-merge(com_738, full_kegg, by.x='gene', by.y='Gene', all.x=T, no.dups=T)
com_862<-merge(com_862, full_kegg, by.x='gene', by.y='Gene', all.x=T, no.dups=T)
com_876<-merge(com_876, full_kegg, by.x='gene', by.y='Gene', all.x=T, no.dups=T)
com_900<-merge(com_900, full_kegg, by.x='gene', by.y='Gene', all.x=T, no.dups=T)
com_908<-merge(com_908, full_kegg, by.x='gene', by.y='Gene', all.x=T, no.dups=T)
com_947<-merge(com_947, full_kegg, by.x='gene', by.y='Gene', all.x=T, no.dups=T)
com_962<-merge(com_962, full_kegg, by.x='gene', by.y='Gene', all.x=T, no.dups=T)

gene_pa2<-merge(gene_pa, full_kegg, by.x='gene', by.y='Gene', all.x=T, no.dups=T)
write.table(gene_pa2, 'full_com.txt', row.names=F, quote=F, sep="\t")

#write to files
#write.table(com_341, "com_341.txt", row.names=F, quote=F, sep="\t")
#write.table(com_738, "com_738-2.txt", row.names=F, quote=F, sep="\t")
#write.table(com_862, "com_862.txt", row.names=F, quote=F, sep="\t")
#write.table(com_876, "com_876.txt", row.names=F, quote=F, sep="\t")
#write.table(com_900, "com_900.txt", row.names=F, quote=F, sep="\t")
#write.table(com_908, "com_908.txt", row.names=F, quote=F, sep="\t")
#write.table(com_947, "com_947.txt", row.names=F, quote=F, sep="\t")
#write.table(com_962, "com_962.txt", row.names=F, quote=F, sep="\t")
#write.table(com_jb, "com_jb.txt", row.names=F, quote=F, sep="\t")



#read back in when added all KO
setwd("/Users/patty/Dropbox/R/cheese_genomes/")
library(readr)
full_kegg <- read_delim("full_kegg.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

full_kegg<-full_kegg[!full_kegg$Level2 == "Endocrine system",]
full_kegg<-full_kegg[!full_kegg$Level2 == "Nervous system",]
full_kegg<-full_kegg[!full_kegg$Level2 == "Nervous system",]
full_kegg<-full_kegg[!full_kegg$Level2 == "Cellular community- eukaryotes",]
full_kegg<-full_kegg[!full_kegg$Level2 == "Cancers: overview",]
full_kegg<-full_kegg[!full_kegg$Level2 == "Cancers: specific types",]
full_kegg<-full_kegg[!full_kegg$Level2 == "Aging",]
full_kegg<-full_kegg[!full_kegg$Level1 == "Diseases",]


full_kegg<-full_kegg[row.names(unique(full_kegg[,c("KO")])),]

library(plyr)
library(tidyr)
com_341<-read.delim("com_341.txt", header=T)
com_341<-merge(com_341, full_kegg, by='KO')
com_341<-com_341 %>% drop_na(Level1)
dim(com_341)

com_738<-read.delim("com_738.txt", header=T)
com_738<-merge(com_738, full_kegg, by='KO')
dim(com_738)

com_862<-read.delim("com_862.txt", header=T)
com_862<-merge(com_862, full_kegg, by='KO')
dim(com_862)

com_876<-read.delim("com_876.txt", header=T)
com_876<-merge(com_876, full_kegg, by='KO')
dim(com_876)

com_962<-read.delim("com_962.txt", header=T)
com_962<-merge(com_962, full_kegg, by='KO')
dim(com_962)

com_908<-read.delim("com_908.txt", header=T)
com_908<-merge(com_908, full_kegg, by='KO')
dim(com_908)

com_900<-read.delim("com_900.txt", header=T)
com_900<-merge(com_900, full_kegg, by='KO')
dim(com_900)

com_947<-read.delim("com_947.txt", header=T)
com_947<-merge(com_947, full_kegg, by='KO')
dim(com_947)

com_jb<-read.delim("com_jb.txt", header=T)
com_jb<-merge(com_jb, full_kegg, by='KO')
com_jb<-com_jb %>% drop_na(Level1)
dim(com_jb)

#add a column for each community name
com_341$com<-rep("Com_341", nrow(com_341))
com_738$com<-rep("Com_738", nrow(com_738))
com_862$com<-rep("Com_862", nrow(com_862))
com_876$com<-rep("Com_876", nrow(com_876))
com_962$com<-rep("Com_962", nrow(com_962))
com_908$com<-rep("Com_908", nrow(com_908))
com_900$com<-rep("Com_900", nrow(com_900))
com_947$com<-rep("Com_947", nrow(com_947))
com_jb$com<-rep("Com_jb", nrow(com_jb))
#merge all together
com_unique<-rbind(com_341, com_738, com_862, com_876, com_962, com_908, com_900, com_947, com_jb)
library(plyr)
com_unique_sum<-ddply(com_unique, c("com", "Level1", "Level2", "Level3"), summarize, Num_genes=length(KO))
library(reshape2)
com_unique_sum_cast<-dcast(com_unique_sum, Level1 + Level2 + Level3 ~ com)
write.table(com_unique_sum_cast, 'KEGG_sum.txt', row.names=F, quote=F, sep='\t')

library(ggplot2)
ggplot(com_unique_sum, aes(com, Level3))+
  geom_tile(aes(fill = Num_genes), colour = "white") + 
  scale_fill_gradient(low = "yellow", high = "green")+
  ylab("KEGG L3")+
  theme_bw()+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(com_unique_sum, aes(com, Level1))+
  geom_tile(aes(fill = Num_genes), colour = "white") + 
  scale_fill_gradient(low = "yellow", high = "green")+
  ylab("KEGG L1")+
  theme_bw()+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(com_unique_sum, aes(com, Level2))+
  geom_tile(aes(fill = Num_genes), colour = "white") + 
  scale_fill_gradient(low = "yellow", high = "green")+
  ylab("KEGG L2")+
  theme_bw()+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
