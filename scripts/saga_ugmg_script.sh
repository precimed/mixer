#!/bin/bash

# Job name:
#SBATCH --job-name=mixer
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

# Job array specification
#SBATCH --array=1-20

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

module load Anaconda3/2019.03
source activate /cluster/home/oleksanf/py3

export MIXER_ROOT=/cluster/projects/nn9114k/oleksanf/github/mixer
export OUTDIR=/cluster/projects/nn9114k/oleksanf/saga/mixer_results/    # must end with a forward slash, /
export SUMSTATnomhc=/cluster/projects/nn9114k/oleksanf/SUMSTAT/TMP/nomhc/
export LDFILE=/cluster/projects/nn9114k/oleksanf/SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld
export BIMFILE=/cluster/projects/nn9114k/oleksanf/SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim
export EXTRACT=/cluster/projects/nn9114k/oleksanf/SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps
export PYTHON=python3
export TRAIT=PGC_SCZ_2014_EUR

$PYTHON $MIXER_ROOT/precimed/mixer.py fit1 \
      --trait1-file $SUMSTATnomhc/$TRAIT.sumstats.gz \
      --out ${OUTDIR}$TRAIT.fit.rep${SLURM_ARRAY_TASK_ID} \
      --bim-file $BIMFILE --ld-file $LDFILE \
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so --threads 20 \
      --extract $EXTRACT \

$PYTHON $MIXER_ROOT/precimed/mixer.py test1 \
      --trait1-file $SUMSTATnomhc/$TRAIT.sumstats.gz \
      --load-params-file $OUT/$TRAIT.fit.rep${SLURM_ARRAY_TASK_ID}.json \
      --out ${OUTDIR}$TRAIT.test.rep${SLURM_ARRAY_TASK_ID} \
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so --threads 20 \
      --bim-file $BIMFILE --ld-file $LDFILE \

