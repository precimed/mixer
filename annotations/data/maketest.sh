#!/bin/bash

CHR_IDS="21 22";
ANNOT_F="../annot.test"
SUMSTATS_F="../sumstats.tets"

for i in ${CHR_IDS}; do
    head chr${i}.bim | cut -f2 > chr${i}.snps.test;
    plink --make-bed --bfile chr${i} --extract chr${i}.snps.test --out ../chr${i}.test;
    plink --bfile ../chr${i}.test --freq --out ../chr${i}.test;
    plink --bfile ../chr${i}.test --r2 gz inter-chr yes-really --ld-window-r2 0.05 --out ../chr${i}.test;
done

cat chr*.snps.test | awk 'BEGIN{OFS="\t";srand(1);print("SNP","ANNOT1","ANNOT2")} {if(rand()<0.5) print($1,1,0); else print($1,0,1);}' > ${ANNOT_F};

echo SNP$'\t'N$'\t'Z > ${SUMSTATS_F} && paste <(cat chr*.snps.test) <(zcat /mnt/seagate10/sumstats/CTG_INTELLIGENCE_2017/CTG_INTELLIGENCE_2017.sumstats.gz | tail -n+2 | shuf | head -n20 | cut -f7,8) | shuf | head -n15 >> ${SUMSTATS_F};

echo Done;

