## Pangenomes
### Convert GBK annotations to GFF annotations for roary
```
mkdir gff_files
#conda install -c bioconda perl-bioperl

#convert with perl script in BioPerl
bp_genbank2gff3.pl -d ./gbk_files -o gff_files
```

## Fix fasta file headers (renamed by fiel name)
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
prokka $i --kingdom Bacteria --cpus 4 --outdir prokka_out_$i --prefix $i
done

#run prokka
for i in SEQ*
do
prokka $i --kingdom Bacteria --cpus 4 --outdir prokka_out_$i --prefix $i
done

#move GFF files
find -type f -name "*.gff" -exec cp {} /home/pattyjk/Desktop/cheese_genomes/gff_files \;
```

## Calculate Pangenome in Roary
```
/home/pattyjk/anaconda3/pkgs/roary-3.7.0-0/bin/roary *.gff -o roary_out
```
