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

for i in *.fna
do
    anvi-gen-contigs-database -f $i -o /home/pattyjk/Desktop/pantoea_ncbi/contigs_db/$i.db -n pantoea
done
