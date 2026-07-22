#!/bin/bash

# fillist <- paste0('/data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb_imp_chr',1:22,'_v3.bgen.bgi')

outpath='/data2/morgante_lab/nklimko/1000G/ukb_snps.txt'

r1='/data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb_imp_chr'
r2='_v3.bgen.bgi'

for n in $(seq 1 22);do
  r3=$r1$n$r2
  # echo $r3
   sqlite3 $r3 'select rsid from variant' >> $outpath
done
