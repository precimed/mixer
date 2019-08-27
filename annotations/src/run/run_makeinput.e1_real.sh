#!/bin/bash
#SBATCH --account=nn9114k
#SBATCH --time=24:00:00 --cpus-per-task 16 --mem-per-cpu=3936M
CFG_FILE="makeinput.e1_real.cfg"

LOG_FILE="cmm_${SLURM_JOBID}.makeinput.log"
CFG_FILE_WITH_ID="cmm_${SLURM_JOBID}.${CFG_FILE}"
#SBATCH --job-name=makeinput_${JOB_ID}
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors
module load plink
module load python3/3.7.0

## Copy files to work directory:
cp -v $SUBMITDIR/../makeinput.py $SCRATCH
cp -v $SUBMITDIR/../cfg/e1_real/${CFG_FILE} $SCRATCH/${CFG_FILE_WITH_ID}
 
## Mark outfiles for automatic copying to $SUBMITDIR:
chkfile ${LOG_FILE} ${CFG_FILE_WITH_ID}

## Run command
cd $SCRATCH
python3 makeinput.py ${CFG_FILE_WITH_ID} > ${LOG_FILE}

