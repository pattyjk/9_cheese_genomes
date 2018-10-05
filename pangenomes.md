## Pangenomes
### Convert GBK annotations to GFF annotations for roary
```
mkdir gff_files
#conda install -c bioconda perl-bioperl

#convert with perl script in BioPerl
bp_genbank2gff3.pl -d ./gbk_files -o gff_files
```

## Calculate Pangenome in Roary
```
/home/pattyjk/anaconda3/pkgs/roary-3.7.0-0/bin/roary *.gff -o roary_out
```
