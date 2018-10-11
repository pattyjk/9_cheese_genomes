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
