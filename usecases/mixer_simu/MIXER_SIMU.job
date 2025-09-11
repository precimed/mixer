#!/bin/bash
#SBATCH --job-name=mixer_simu
#SBATCH --output=MIXER_SIMU-%A_%a.txt
#SBATCH --error=MIXER_SIMU-%A_%a.txt
#SBATCH --open-mode=truncate
#SBATCH --account=p697_norment
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8000M
#SBATCH --array=1-20

module purge
module load singularity/3.7.1

export COMORMENT=/ess/p697/data/durable/s3-api/github/comorment
export SINGULARITY_BIND="$COMORMENT/mixer/reference:/REF:ro"
export SIF=$COMORMENT/mixer/singularity
export MIXER_COMMON_ARGS="--chr2use 21 --ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim"
export REP="rep${SLURM_ARRAY_TASK_ID}"
export EXTRACT="--extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.$REP.snps"

export PYTHON="singularity exec --home $PWD:/home $SIF/mixer.sif python"

$PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file unique.1.sumstats --out unique.1.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file unique.2.sumstats --out unique.2.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py fit2 $MIXER_COMMON_ARGS $EXTRACT --trait1-file unique.1.sumstats --trait2-file unique.2.sumstats --trait1-params unique.1.fit.$REP.json --trait2-params unique.2.fit.$REP.json --out unique.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py test2 $MIXER_COMMON_ARGS --trait1-file unique.1.sumstats --trait2-file unique.2.sumstats --load-params unique.fit.$REP.json --out unique.test.$REP

$PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file shared.1.sumstats --out shared.1.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file shared.2.sumstats --out shared.2.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py fit2 $MIXER_COMMON_ARGS $EXTRACT --trait1-file shared.1.sumstats --trait2-file shared.2.sumstats --trait1-params shared.1.fit.$REP.json --trait2-params shared.2.fit.$REP.json --out shared.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py test2 $MIXER_COMMON_ARGS --trait1-file shared.1.sumstats --trait2-file shared.2.sumstats --load-params shared.fit.$REP.json --out shared.test.$REP

$PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file partial.1.sumstats --out partial.1.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file partial.2.sumstats --out partial.2.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py fit2 $MIXER_COMMON_ARGS $EXTRACT --trait1-file partial.1.sumstats --trait2-file partial.2.sumstats --trait1-params partial.1.fit.$REP.json --trait2-params partial.2.fit.$REP.json --out partial.fit.$REP
$PYTHON /tools/mixer/precimed/mixer.py test2 $MIXER_COMMON_ARGS --trait1-file partial.1.sumstats --trait2-file partial.2.sumstats --load-params partial.fit.$REP.json --out partial.test.$REP
