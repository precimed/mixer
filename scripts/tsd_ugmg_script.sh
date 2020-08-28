#!/bin/bash

# Job name:
#SBATCH --job-name=mixer
#
# Project:
#SBATCH --account=p33_norment
#
# Wall clock limit:
#SBATCH --time=16:00:00
#
#SBATCH --cpus-per-task=8

# Max memory usage:
#SBATCH --mem-per-cpu=7600M

# Job array specification
#SBATCH --array=1-20

## Set up job environment:
source /cluster/bin/jobsetup
set -o errexit

#module init
module load CMake/3.15.3-GCCcore-8.3.0
module load Boost/1.73.0-GCCcore-8.3.0
module load Python/3.7.4-GCCcore-8.3.0
source /cluster/projects/p33/users/ofrei/py3/bin/activate

export MIXER_ROOT=/cluster/projects/p33/users/ofrei/github/mixer
export OUTDIR=/cluster/projects/p33/users/ofrei/mixer_results/    # must end with a forward slash, /
export SUMSTATnomhc=/cluster/projects/p33/users/ofrei/SUMSTAT/TMP/nomhc/
export LDFILE=/cluster/projects/p33/users/ofrei/SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld
export BIMFILE=/cluster/projects/p33/users/ofrei/SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim
export EXTRACT=/cluster/projects/p33/users/ofrei/SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps
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

