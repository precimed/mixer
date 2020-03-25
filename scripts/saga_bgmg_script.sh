#!/bin/bash

# **************************************************************
# STEP 1: ABEL specifications for job submission
# **************************************************************

# Job name:
#SBATCH --job-name=ofrebgmg
#
# Project:
#SBATCH --account=nn9114k
#
# Wall clock limit:
#SBATCH --time=1-00:00:00
#
#SBATCH --cpus-per-task=20

# Max memory usage:
#SBATCH --mem-per-cpu=4600M

# Email about job status
#SBATCH --mail-type=ALL

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

module load Anaconda3/2019.03
source activate /cluster/home/oleksanf/py3

export MIXER_ROOT=/cluster/projects/nn9114k/oleksanf/github/mixer_plsa
export SUMSTAT=/cluster/projects/nn9114k/oleksanf/SUMSTAT
export OUT_DIR=/cluster/projects/nn9114k/oleksanf/saga/mixer_test/saga2
export PYTHON=python3
export TRAIT1=PGC_SCZ_2014_EUR
export TRAIT2=SSGAC_EDU_2018_no23andMe

$PYTHON $MIXER_ROOT/precimed/mixer.py fit \
      --trait1-file $SUMSTAT/TMP/nomhc/$TRAIT1.sumstats.gz \
      --trait2-file $SUMSTAT/TMP/nomhc/$TRAIT2.sumstats.gz \
      --trait1-params-file $OUT_DIR/$TRAIT1.fit.json \
      --trait2-params-file $OUT_DIR/$TRAIT2.fit.json \
      --out $OUT_DIR/${TRAIT1}_vs_${TRAIT2}.fit \
      --extract $SUMSTAT/LDSR/w_hm3.justrs --ci-alpha 0.05 \
      --bim-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --frq-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.frq \
      --plink-ld-bin0 $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin \
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so \

#      --chr2use 1 --diffevo-fast-repeats 2 --fit-sequence diffevo-fast neldermead-fast

$PYTHON $MIXER_ROOT/precimed/mixer.py fit \
      --trait1-file $SUMSTAT/TMP/nomhc/$TRAIT1.sumstats.gz \
      --trait2-file $SUMSTAT/TMP/nomhc/$TRAIT2.sumstats.gz \
      --load-params-file $OUT_DIR/${TRAIT1}_vs_${TRAIT2}.fit.json \
      --out $OUT_DIR/${TRAIT1}_vs_${TRAIT2}.test \
      --bim-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --frq-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.frq \
      --plink-ld-bin0 $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin \
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so \
      --fit-sequence load inflation --qq-plots --kmax 100 \
 
#     --chr2use 1

