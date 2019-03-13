## Pangenomes with Roary
### Convert GBK annotations to GFF annotations for Roary
```
mkdir gff_files
#conda install -c bioconda perl-bioperl

#convert with perl script in BioPerl
bp_genbank2gff3.pl -d ./gbk_files -o gff_files
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

#run prokka
for i in SEQ*
do
prokka $i --kingdom Bacteria --cpus 4 --outdir prokka_out_$i --prefix $i --locustag $i
done

#move GFF files
find -type f -name "*.gff" -exec cp {} /home/pattyjk/Desktop/cheese_genomes/gff_files \;
```

## Calculate Pangenome in Roary
```
/home/pattyjk/anaconda3/pkgs/roary-3.7.0-0/bin/roary *.gff -o roary_out
```

## Extract gene names
```
for ID in cat $HOME/Dropbox/R/cheese_genomes/unique_gene_list.txt; do grep $ID com_341.fna.gff; done > genes_341.txt

for ID in cat $HOME/Dropbox/R/cheese_genomes/unique_gene_list.txt; do grep $ID com_862.fna.gff; done > genes_862.txt

for ID in cat $HOME/Dropbox/R/cheese_genomes/unique_gene_list.txt; do grep $ID com_900.fna.gff; done > genes_900.txt

for ID in cat $HOME/Dropbox/R/cheese_genomes/unique_gene_list.txt; do grep $ID com_947.fna.gff; done > genes_947.txt

for ID in cat $HOME/Dropbox/R/cheese_genomes/unique_gene_list.txt; do grep $ID com_JB.fna.gff; done > genes_jb.txt

for ID in cat $HOME/Dropbox/R/cheese_genomes/unique_gene_list.txt; do grep $ID com_738.fna.gff; done > genes_738.txt

for ID in cat $HOME/Dropbox/R/cheese_genomes/unique_gene_list.txt; do grep $ID com_876.fna.gff; done > genes_876.txt

for ID in cat $HOME/Dropbox/R/cheese_genomes/unique_gene_list.txt; do grep $ID com_908.fna.gff; done > genes_908.txt

for ID in cat $HOME/Dropbox/R/cheese_genomes/unique_gene_list.txt; do grep $ID com_962.fna.gff; done > genes_962.txt

mkdir genes
mv genes_* genes
cd genes
cat genes_* > all_genes.txt
```

## Get gene names
```
#get database (all bact/arch) from NCBI ftp
ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Archaea_Bacteria/

#extract all names
for ID in cat /home/pattyjk/Dropbox/R/cheese_genomes/all_genes.txt; do grep $ID /home/pattyjk/Desktop/All_Archaea_Bacteria.gene_info; done > Desktop/gene_info.txt 
```
