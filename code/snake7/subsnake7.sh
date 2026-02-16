#!/bin/bash
savePath='/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_nklimko/11_snake/'

function queuePrint {
   sbatch '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/snake7/snake7.sbatch' | grep -o -E '[0-9]+'
#echo temp
}

ticket=$(queuePrint)
savePath+=$ticket

export NN7=$savePath
