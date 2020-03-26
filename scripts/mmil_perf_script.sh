
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

#$ -l h_vmem=120G
#$ -l h_rt=36:00:00
##$ -pe dmp4 16

##$ -l h=!(mmil-compute-5-7.local|mmil-compute-5-14.local|mmil-compute-5-5.local|mmil-compute-5-2.local|mmil-compute-5-10.local|mmil-compute-8-1|mmil-compute-8-2|mmil-compute-8-3|mmil-compute-8-4|mmil-compute-7-0|mmil-compute-7-1)

set MIXER_ROOT=/home/oleksandr/github/mixer
set SUMSTAT=/space/syn03/1/data/GWAS/SUMSTAT
set OUT_DIR=/space/syn03/1/data/oleksandr/mixer_test/perf_cluster7
set PYTHON=/home/oleksandr/miniconda3/bin/python3
set TRAIT1=PGC_SCZ_2014_EUR
set TRAIT2=SSGAC_EDU_2018_no23andMe

##$PYTHON $MIXER_ROOT/precimed/mixer.py perf \
##      --trait1-file $SUMSTAT/TMP/ldsr/$TRAIT1.sumstats.gz \
##      --out $OUT_DIR/${TRAIT1} \
##      --bim-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
##      --plink-ld-bin $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
##      --lib $MIXER_ROOT/src/build/lib/libbgmg.so \
##      --kmax 20000 2000 200 \
##      --threads 1 2 4 8 12 16 20  30 40  \
##
##      --chr2use 21 

$PYTHON $MIXER_ROOT/precimed/mixer.py perf \
      --trait1-file $SUMSTAT/TMP/ldsr/$TRAIT1.sumstats.gz \
      --trait2-file $SUMSTAT/TMP/ldsr/$TRAIT2.sumstats.gz \
      --out $OUT_DIR/${TRAIT1}_vs_${TRAIT2} \
      --bim-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --plink-ld-bin $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so \
      --kmax 20000 2000 200 \
      --threads 1 2 4 8 12 16 20  30 40  \

