## Pangenomes with Roary
### Convert GBK annotations to GFF annotations for Roary
```
mkdir gff_files
#conda install -c bioconda perl-bioperl

#convert with perl script in BioPerl
bp_genbank2gff3.pl *.gbk --outdir gf_files

#run Roary
/home/pattyjk/anaconda3/pkgs/roary-3.7.0-0/bin/roary *.gff -o roary_out
```

## Fix fasta file headers for Prokka (renamed by file name)
```
for file in *.fa;
do
sed -i "s/>.*/>${file%%.*}/" "$file" ;
done

perl -i.bak -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' *.fa
perl -i.bak -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' SEQ*
rm *.bak
```

## Prokka
```
#run prokka
for i in *.fa 
do
prokka $i --kingdom Bacteria --cpus 4 --outdir prokka_out_$i --prefix $i --locustag $i
done

mkdir gff_files

#move GFF files
find -type f -name "*.gff" -exec cp {} /home/pattyjk/9_cheese_genomes/assembled_genomes/gff_files \;
```

## Calculate Pangenome in Roary
```
cd gff_files
/home/pattyjk/anaconda3/pkgs/roary-3.7.0-0/bin/roary *.gff -o roary_out
```

## Analyze pangenome in R
```
library(plyr)
library(ggplot2)
library(vegan)
library(reshape2)
library(dplyr)
library(tidyr)

#read in data
pan_full<-read.delim("9_cheese_genomes/roary_out/gene_presence_absence.csv", sep=',', header=T)

#grab columns of interest
pan_full2<-pan_full[,c(1,3, 15:41)]

#create coumn with the presence in all genomes
pan_full2<-unite(pan_full2, "full_genomes", c(names(pan_full2[,3:29])), sep=',', remove=T)

#remove excess commas
pan_full2$full_genomes<-gsub(",,,", "", pan_full2$full_genomes)


#split by community
com_JB<-select(pan_full, matches('_JB'))
com_341<-select(pan_full, matches('_341_'))
com_738<-select(pan_full, matches('_738_'))
com_862<-select(pan_full, matches('_862_'))
com_876<-select(pan_full, matches('_876_'))
com_900<-select(pan_full, matches('_900_'))
com_908<-select(pan_full, matches('_908_'))
com_947<-select(pan_full, matches('_947_'))
com_962<-select(pan_full, matches('_962_'))

#create column for each community that has the presence of each gene
com_JB<-unite(com_JB, 'com_JB_genomes', c(names(com_JB)), sep=',', remove=F)
com_341<-unite(com_341, 'com_341_genomes', c(names(com_341)), sep=',', remove=F)
com_738<-unite(com_738, 'com_738_genomes', c(names(com_738)), sep=',', remove=F)
com_862<-unite(com_862, 'com_862_genomes', c(names(com_862)), sep=',', remove=F)
com_876<-unite(com_876, 'com_876_genomes', c(names(com_876)), sep=',', remove=F)
com_900<-unite(com_900, 'com_900_genomes', c(names(com_900)), sep=',', remove=F)
com_908<-unite(com_908, 'com_908_genomes', c(names(com_908)), sep=',', remove=F)
com_947<-unite(com_947, 'com_947_genomes', c(names(com_947)), sep=',', remove=F)
com_962<-unite(com_962, 'com_962_genomes', c(names(com_962)), sep=',', remove=F)

#catenate gene communities
com_full<-cbind(com_JB, com_341, com_738, com_862, com_876, com_900, com_908, com_947, com_962)

#bring data back together
pan_full3<-cbind(pan_full, select(com_full, matches('_genomes')))


#read in presence/absence data
pan<-read.delim("9_cheese_genomes/roary_out/gene_presence_absence.Rtab", header=T)

#reshape and fix genome names
pan_m<-melt(pan)
pan_m$variable<-gsub("Brachy_", "", pan_m$variable)
pan_m$variable<-gsub("SEQ_", "", pan_m$variable)
pan_m$variable<-gsub("Brevi_", "", pan_m$variable)
pan_m$variable<-gsub("_[[:digit:]]+.fa", "", pan_m$variable)
pan_m$variable<-gsub(".fa", "", pan_m$variable)
pan_m$variable<-gsub("JB7", "JB5", pan_m$variable)
pan_m$variable<-gsub("BC9", "JB5", pan_m$variable)

#summarize data
pan_sum<-ddply(pan_m, c("variable", 'Gene'), summarize, total=sum(value))

ggplot(pan_sum, aes(Gene, variable, fill=total))+
geom_tile()

#recast data
pan_cast<-dcast(pan_sum, Gene ~ variable, value.var = "total")

#create presence absence data frame for genes
pan_cast_pa<-pan_cast[,-1]
pan_cast_pa<-decostand(pan_cast_pa, method='pa')
names(pan_cast_pa)<-c("pa_341", "pa_738", "pa_862", "pa_876", "pa_900", "pa_908", "pa_947", "pa_962", "pa_JB5")

#bind data back together
pan_cast2<-cbind(pan_cast, pan_cast_pa)

#merge all data back together
pangenomes<-merge(pan_cast2, pan_full3, by='Gene')

#write data to file
write.table(pangenomes, 'pangenomes.txt', sep='\t', row.names=F, quote=F)
```

## Identify genes unique to each community
```
#make a table of unique genes
panny<-pangenomes[,11:19]
unique_gene_list<-which(rowSums(panny) == 1)
names(unique_gene_list)<-"RowSum"
dim(unique_gene_list)
#6250 genes

#extract unique genes from master table
unique_genes<-subset(pangenomes, row.names(pangenomes) %in% unique_gene_list)
dim(unique_genes)
#6250 by 68

#create tables for each community
com_341_unique<- unique_genes[unique_genes$pa_341 == 1,] %>% select("341", "pa_341", "Gene", 'Annotation', 'com_341_genomes')
nrow(com_341_unique)
#569

com_738_unique<- unique_genes[unique_genes$pa_738 == 1,] %>% select("738", "pa_738", "Gene", 'Annotation', 'com_738_genomes')
nrow(com_738_unique)
#1014

com_862_unique<- unique_genes[unique_genes$pa_862 == 1,] %>% select("862", "pa_862", "Gene", 'Annotation', 'com_862_genomes')
nrow(com_862_unique)
#702

com_876_unique<- unique_genes[unique_genes$pa_876 == 1,] %>% select("876", "pa_876", "Gene", 'Annotation', 'com_876_genomes')
nrow(com_876_unique)
#496

com_900_unique<- unique_genes[unique_genes$pa_900 == 1,] %>% select("900", "pa_900", "Gene", 'Annotation', 'com_900_genomes')
nrow(com_900_unique)
#774

com_908_unique<- unique_genes[unique_genes$pa_908 == 1,] %>% select("908", "pa_908", "Gene", 'Annotation', 'com_908_genomes')
nrow(com_908_unique)
#484

com_947_unique<- unique_genes[unique_genes$pa_947 == 1,] %>% select("947", "pa_947", "Gene", 'Annotation', 'com_947_genomes')
nrow(com_947_unique)
#577

com_962_unique<- unique_genes[unique_genes$pa_962 == 1,] %>% select("962", "pa_962", "Gene", 'Annotation', 'com_962_genomes')
nrow(com_962_unique)
#407

com_JB5_unique<- unique_genes[unique_genes$pa_JB5 == 1,] %>% select("JB5", "pa_JB5", "Gene", 'Annotation', 'com_JB_genomes')
nrow(com_JB5_unique)
#1227
```
