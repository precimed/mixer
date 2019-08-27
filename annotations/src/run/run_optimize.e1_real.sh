#!/bin/bash
#SBATCH --account=nn9114k
#SBATCH --time=24:00:00 --cpus-per-task 16 --mem-per-cpu=3936M
trait_id=$1
i_repeat=$2
CFG_FILE="optimize.e1_real.template.${trait_id}.${i_repeat}.cfg"

LOG_FILE="cmm_${SLURM_JOBID}.optimize.log"
CFG_FILE_WITH_ID="cmm_${SLURM_JOBID}.${CFG_FILE}"
#SBATCH --job-name=optimize_${JOB_ID}
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors
module load python3/3.7.0
module load gsl/2.2

## Copy files to work directory:
mkdir $SCRATCH/c
cp -v $SUBMITDIR/../optimize.py $SCRATCH
cp -v $SUBMITDIR/../c/__init__.py $SCRATCH/c
cp -v $SUBMITDIR/../c/_cmmcost_omp.c $SCRATCH/c
cp -v $SUBMITDIR/../c/_cmmcost_omp.h $SCRATCH/c
cp -v $SUBMITDIR/../c/setup_cmmcost_omp_abel.py $SCRATCH/c
cp -v $SUBMITDIR/../c/cmmcost_omp.pyx $SCRATCH/c
cp -v $SUBMITDIR/../cfg/e1_real/${CFG_FILE} $SCRATCH/${CFG_FILE_WITH_ID}
 
## Mark outfiles for automatic copying to $SUBMITDIR:
chkfile ${LOG_FILE} ${CFG_FILE_WITH_ID}

## Compile
cd $SCRATCH/c
python3 setup_cmmcost_omp_abel.py build_ext --inplace

## Run command
cd $SCRATCH
python3 optimize.py ${CFG_FILE_WITH_ID} > ${LOG_FILE}

