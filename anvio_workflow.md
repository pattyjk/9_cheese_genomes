## Anvio workflow
```
#load anvio with Anaconda
source activate anvio5

#fix fasta files to be compatible with Anvio
mkdir fixed_reps
mkdir fixed_fasta

#run for loop
for i in *
do
    anvi-script-reformat-fasta $i -o fixed_fasta/$i --simplify-names --report-file=fixed_reps/$i
done

#make contigs database for each genome
mkdir contigs_db
cd fixed_fasta

for i in *
do
    anvi-gen-contigs-database -f $i -o /home/pattyjk/Desktop/cheese_genomes/anvio_work/$i.db -n pantoea
done

#generate COGs
cd ..
cd contigs_db

#make a contigis database if not already done (anvi-setup-ncbi-cogs --num-threads 6)
for i in *.db
do
anvi-run-ncbi-cogs -c $i --num-threads 8
done

#run HMMs for single copy genes
for i in *.db
do
anvi-run-hmms -c $i -T 8
done

cd ..
```

## Make database
```
ls -R1 ./home/pattyjk/Desktop/cheese_genomes/anvio_work/contigs_db | while read l; do case $l in *:) d=${l%:};; "") d=;; *) echo "$d/$l";; esac; done > /home/pattyjk/Desktop/cheese_genomes/anvio_work/db_files.txt

cd /home/pattyjk/Desktop/cheese_genomes/anvio_work

sed 's/\///g' db_files.txt > gen_names.txt

anvi_gen<-read.delim('db_files.txt', header=F)
names(anvi_gen)<-'contigs_db_path'

#read in genome names
anvi_names<-read.delim('gen_names.txt', header=F)
names(anvi_names)<-'name'

#catenate
anvi_gen<-cbind(anvi_names, anvi_gen)

#write file
write.table(anvi_gen, 'anvi_gen.txt', quote=F, row.names=F, sep='\t')

#exit R
quit()
n
```

## Catenate genomes in a DB file
```
anvi-gen-genomes-storage -e anvi_gen.txt -o cheese-GENOMES.db 

#pangenome analysis
anvi-pan-genome -g cheese-GENOMES.db -n cheese --enforce-hierarchical-clustering --min-occurrence 2

#view analysis
anvi-display-pan -g cheese-GENOMES.db -p cheese/cheese-PAN.db
```

## Phylogenomics on RNA proteins
```
#generate genome tree
anvi-get-sequences-for-hmm-hits --external-genomes anvi_gen.txt -o concatenated-proteins.fa --hmm-source Campbell_et_al --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 --return-best-hit --get-aa-sequence --concatenate

anvi-gen-phylogenomic-tree -f concatenated-proteins.fa -o pantoea_tree.tree
