## Example of GSA-MiXeR analysis

See [scripts/GSA_MIXER.job](../scripts/GSA_MIXER.job) for an example of script running GSA-MiXeR analysis.
Optionally, re-format the results using [process_gsa_mixer_output.py](scripts/process_gsa_mixer_output.py) script.

### Slurm configuration

GSA-MiXeR analysis takes around 6 to 12 hours using 8-core machine, and utilize around 30 GB of RAM. The following configuration might be a good starting point:
```
#SBATCH --job-name=gsamixer
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account=nn9114k
#SBATCH --mem-per-cpu=4000M
#SBATCH --cpus-per-task=8     # if you change keep --threads $THREADS argument in sync
```

### Define data location, and singularity-related variables

Check [scripts/GSA_MIXER.job](../scripts/GSA_MIXER.job) for more up-to-date version.
```
export GITHUB=/cluster/projects/nn9114k/github
export MIXER_SIF=${GITHUB}/precimed/gsa-mixer/containers/latest/mixer.sif
export SUMSTATS_FOLDER=/cluster/projects/nn9114k/oleksanf/gsa-mixer/sumstats
export SUMSTATS_FILE=PGC_SCZ_0518_EUR
export OUT_FOLDER=/cluster/projects/nn9114k/oleksanf/gsa-mixer/out2
export BIND="--bind /cluster/projects/nn9114k:/cluster/projects/nn9114k"

export REFERENCE_FOLDER=${GITHUB}/comorment/mixer/reference
export BIM_FILE=${GITHUB}/comorment/mixer/reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim
export LOADLIB_FILE=${GITHUB}/comorment/mixer/reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin
export ANNOT_FILE=${GITHUB}/comorment/mixer/reference/ldsc/1000G_EUR_Phase3_plink/baseline_v2.2_1000G.EUR.QC.@.annot.gz

export MAGMA_BFILE=${GITHUB}/comorment/mixer/reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@
export MAGMA_GENE_LOC=${GITHUB}/comorment/mixer/reference/magma-gene-annot_10mar2023.csv
export MAGMA_SET_ANNOT=${GITHUB}/comorment/mixer/reference/magma-geneset-annot_10mar2023.csv

export PYTHON="singularity exec --home pwd:/home $BIND ${MIXER_SIF} python"
export MIXER_PY="$PYTHON /tools/mixer/precimed/mixer.py"
export MAGMA="singularity exec --home pwd:/home $BIND ${MIXER_SIF} magma"
```

### Perform GSA-MIXER analysis

Define a few configuration options shared between baseline and full models:
```
export EXTRA_FLAGS="--seed 1000 --exclude-ranges MHC --hardprune-r2 0.6 --threads 8 "
```
It's recommended to update the ``--threads`` argument so that it is in sync with SLURM's ``--cpus-per-task``.

Baseline model:
```
${MIXER_PY} plsa --gsa-base \
        --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.chr@.sumstats.gz \
        --out ${OUT_FOLDER}/${SUMSTATS_FILE}_base \
        --bim-file ${BIM_FILE} --use-complete-tag-indices --loadlib-file ${LOADLIB_FILE} \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_10mar2023.csv \
        --annot-file ${ANNOT_FILE} \
        ${EXTRA_FLAGS}
```

Enrichment model:
```
${MIXER_PY} plsa --gsa-full \
        --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.chr@.sumstats.gz \
        --out ${OUT_FOLDER}/${SUMSTATS_FILE}_full \
        --bim-file ${BIM_FILE} --use-complete-tag-indices --loadlib-file ${LOADLIB_FILE} \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-gene-annot_10mar2023.csv \
        --go-file-test ${REFERENCE_FOLDER}/gsa-mixer-hybridLOO-annot_10mar2023.csv \
        --annot-file ${ANNOT_FILE} \
        --load-params-file ${OUT_FOLDER}/${SUMSTATS_FILE}_base.json \
        ${EXTRA_FLAGS}
```

Key output file is ``${OUT_FOLDER}/${SUMSTATS_FILE}_full.go_test_enrich.csv``.
You may check ``out_example`` folder for a pre-generated example of such file.

### Perform MAGMA analysis

```
# MAGMA analysis - annotate snps to genes
$MAGMA --snp-loc ${MAGMA_BFILE}.bim \
       --gene-loc ${MAGMA_GENE_LOC} \
       --out ${OUT_FOLDER}/${SUMSTATS_FILE}_magma.step1 \
       --annotate window=10

# MAGMA analysis - compute gene-level p-values
$MAGMA --pval ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats snp-id=RSID pval=P ncol=N \
       --bfile ${MAGMA_BFILE} \
       --gene-annot ${OUT_FOLDER}/${SUMSTATS_FILE}_magma.step1.genes.annot \
       --out ${OUT_FOLDER}/${SUMSTATS_FILE}_magma.step2

# MAGMA analysis - compute geneset-level p-values
$MAGMA --gene-results ${OUT_FOLDER}/${SUMSTATS_FILE}_magma.step2.genes.raw \
       --set-annot ${MAGMA_SET_ANNOT} \
       --out ${OUT_FOLDER}/${SUMSTATS_FILE}_magma
```

Key output files are
``${OUT_FOLDER}/${SUMSTATS_FILE}_magma.step2.genes.out`` and ``${OUT_FOLDER}/${SUMSTATS_FILE}_magma.gsa.out``.
You may check ``out_example`` folder for a pre-generated example of such files.

### Re-format the results

[process_gsa_mixer_output.py](scripts/process_gsa_mixer_output.py) script can be used to re-format the results,
assuming that GSA-MiXeR and MAGMA outputs are stored with ``${OUT_FOLDER}/${SUMSTATS_FILE}_full`` and ``${OUT_FOLDER}/${SUMSTATS_FILE}_magma`` prefixes, respectively.
We include a few sample files in this repository allowing to test out the [scripts/process_gsa_mixer_output.py](scripts/process_gsa_mixer_output.py) script.
You will need to make a trivial changes to the script, pointing it to the location of the input files.
To run the script you may use python installation from MiXeR's docker or singularity container, as it includes all required dependencies:
```
export PYTHON="singularity exec --home pwd:/home $BIND ${MIXER_SIF} python"
$PYTHON process_gsa_mixer_output.py
```
This yields ``SupplementaryTables.xlsx`` file, with columns named as in respective supplementary tables from the GSA-MiXeR's publication.
