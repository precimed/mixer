#!/bin/bash

# Execute job in the queue "std.q" unless you have special requirements.
##$ -q std.q
##$ -q all_24.q

# The SGE batch system uses the current directory as working directory.
# Both files (output.dat and error.dat) will be placed in the current
# directory. The batch system assumes to find the executable in this directory.
#$ -cwd

# Redirect output stream to this file.
##$ -o output.dat

# Redirect error stream to this file.
##$ -e error.dat

# Send status information to this email address.
##$ -M oleksandr.frei@gmail.com

# Send an e-mail when the job is done.
##$ -m e

# Use bash
##$ -S /bin/bash

#$ -l h_vmem=120G
#$ -l h_rt=36:00:00
##$ -pe dmp4 16

##$ -l h=!(mmil-compute-5-7.local|mmil-compute-5-14.local|mmil-compute-5-5.local|mmil-compute-5-2.local|mmil-compute-5-10.local|mmil-compute-8-1|mmil-compute-8-2|mmil-compute-8-3|mmil-compute-8-4|mmil-compute-7-0|mmil-compute-7-1)

## https://unix.stackexchange.com/questions/277981/gnu-parallel-immediately-display-job-stderr-stdout-one-at-a-time-by-jobs-order

set MIXER_ROOT=/home/oleksandr/github/mixer 
set SUMSTAT=/space/syn03/1/data/GWAS/SUMSTAT
set OUT_DIR=/space/syn03/1/data/oleksandr/mixer_test/cleanup_legacy
set PYTHON=/home/oleksandr/miniconda3/bin/python3
set TRAIT=PGC_SCZ_2014_EUR
#set TRAIT=SSGAC_EDU_2018_no23andMe 

$PYTHON $MIXER_ROOT/precimed/mixer.py fit \
      --trait1-file $SUMSTAT/TMP/ldsr/$TRAIT.sumstats.gz \
      --out $OUT_DIR/$TRAIT.fit \
      --extract $SUMSTAT/LDSR/w_hm3.justrs --ci-alpha 0.05 \
      --bim-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --frq-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.frq \
      --plink-ld-bin0 $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin \
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so \

#      --chr2use 1 --diffevo-fast-repeats 2 --fit-sequence diffevo-fast neldermead-fast

$PYTHON $MIXER_ROOT/precimed/mixer.py fit \
      --trait1-file $SUMSTAT/TMP/ldsr/$TRAIT.sumstats.gz \
      --load-params-file $OUT_DIR/$TRAIT.fit.json \
      --out $OUT_DIR/$TRAIT.test \
      --bim-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --frq-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.frq \
      --plink-ld-bin0 $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin \
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so \
      --fit-sequence load inflation --power-curve --qq-plots --kmax 100 \

#     --chr2use 1
