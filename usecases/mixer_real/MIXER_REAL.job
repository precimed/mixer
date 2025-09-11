#!/bin/bash
#SBATCH --job-name=mixer_real
#SBATCH --output=MIXER_REAL-%A_%a.txt
#SBATCH --error=MIXER_REAL-%A_%a.txt
#SBATCH --open-mode=truncate
#SBATCH --account=p697_norment
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8000M
#SBATCH --array=1-20

module purge
module load singularity/3.7.1

export COMORMENT=/ess/p697/data/durable/s3-api/github/comorment
export SINGULARITY_BIND="$COMORMENT/mixer/reference:/REF:ro"
export SIF=$COMORMENT/mixer/singularity
export MIXER_COMMON_ARGS="--ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --threads 16"
export REP="rep${SLURM_ARRAY_TASK_ID}"
export EXTRACT="--extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.$REP.snps"

export PYTHON="singularity exec --home=$PWD:/home $SIF/mixer.sif python"

$PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file SCZ.sumstats.gz --out SCZ.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file INT.sumstats.gz --out INT.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py fit2 $MIXER_COMMON_ARGS $EXTRACT --trait1-file SCZ.sumstats.gz --trait2-file INT.sumstats.gz --trait1-params SCZ.fit.$REP.json --trait2-params INT.fit.$REP.json --out SCZ_vs_INT.fit.$REP

$PYTHON /tools/mixer/precimed/mixer.py test1 $MIXER_COMMON_ARGS --trait1-file SCZ.sumstats.gz --load-params SCZ.fit.$REP.json --out SCZ.test.$REP
$PYTHON /tools/mixer/precimed/mixer.py test1 $MIXER_COMMON_ARGS --trait1-file INT.sumstats.gz --load-params INT.fit.$REP.json --out INT.test.$REP
$PYTHON /tools/mixer/precimed/mixer.py test2 $MIXER_COMMON_ARGS --trait1-file SCZ.sumstats.gz --trait2-file INT.sumstats.gz --load-params SCZ_vs_INT.fit.$REP.json --out SCZ_vs_INT.test.$REP
