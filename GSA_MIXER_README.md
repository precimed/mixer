## Introduction

GSA-MiXeR is a new technique for competitive gene-set analysis, which fits a model for gene-set heritability enrichments for complex human traits, thus allowing the quantification of partitioned heritability and fold enrichment for small gene-sets.

Cross-trait MiXeR is a statistical tool which quantifies polygenic overlap between complex traits irrespective of genetic correlation, using GWAS summary statistics. MiXeR results are presented as a Venn diagram of unique and shared polygenic components across traits.

For a real-world application of the GSA-MiXeR you will need to perform the following steps:
* [Install GSA-MiXeR](#install-gsa-mixer) using Docker or singularity containers
* (optionally) test our your installation using [Getting Started Example](#getting-started-example) with tiny dummy data
* Format your GWAS summary statistics according to [Input Data Format](#input-data-formats)
* [Download pre-generated LD matrix and other reference files](#download-ld-matrix-and-other-reference-files) based on 1000 Genomes European population
* [Perform GSA-MiXeR and MAGMA analyses](#perform-gsa-mixer-and-magma-analyses) using modified version of the [GSA_MIXER.job](scripts/GSA_MIXER.job) script; optionally, re-format the results using [process_gsa_mixer_output.py](scripts/process_gsa_mixer_output.py) script.

GSA-MiXeR software also supports univariate and bivariate (cross-trait) MiXeR analyses from [Frei et al, 2019](https://www.nature.com/articles/s41467-019-10310-0). This replaces previous ``MiXeR v1.3`` software package (https://github.com/precimed/mixer) and its pre-built singularity container (https://github.com/comorment/mixer). The [scripts/MIXER.job](scripts/MIXER.job) script give an example of performing univariate and bivariate analyses using GSA-MiXeR software, inline with procedure previously employed by ``MiXeR v1.3``.

For further information refer to [Command-line reference](#command-line-reference) section.
We also provide instructions on how to [generate your own LD reference](#generate-ld-reference) files, for example using UKB or HRC genotypes.

The history of software changes is available in the [CHANGELOG.md](CHANGELOG.md) file.

Please cite relevant publications if you use MiXeR software in your research work.

* for univariate analysis: D. Holland et al., Beyond SNP Heritability: Polygenicity and Discoverability Estimated for Multiple Phenotypes with a Univariate Gaussian Mixture Model, PLOS Genetics, 2020, https://doi.org/10.1371/journal.pgen.1008612
* for cross-trait analysis: O.Frei et al., Bivariate causal mixture model quantifies polygenic overlap between complex traits beyond genetic correlation, Nature Communications, 2019, https://www.nature.com/articles/s41467-019-10310-0
* for gene-set enrichment analysis: O.Frei et al., Improved functional mapping with GSA-MiXeR implicates biologically specific gene-sets and estimates enrichment magnitude, Nature Genetics, 2024, https://www.nature.com/articles/s41588-024-01771-1

## Install GSA-MiXeR

GSA-MiXeR software is released both as a Docker container, and as a pre-compiled singularity (apptainer) container. Use the following commands to check if you have Docker and/or singulrity available in your environment:
```
# check if Docker software is installed
>docker --version
Docker version 20.10.7, build 20.10.7-0ubuntu5~21.04.2

# check if singularity software is installed
>singularity --version
singularity version 3.7.4
```

To dowload Docker version of the GSA-MiXeR, use the following command:
```
docker pull ghcr.io/precimed/gsa-mixer:latest
export DOCKER_RUN="docker run -v $PWD:/home -w /home"
export MIXER_PY="$DOCKER_RUN ghcr.io/precimed/gsa-mixer:latest python /tools/mixer/precimed/mixer.py"
```
To download singularity version of the GSA-MiXeR, use the following command:
```
oras pull ghcr.io/precimed/gsa-mixer_sif:latest
export MIXER_SIF=<path>/gsa-mixer.sif
export MIXER_PY="singularity exec --home pwd:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"
```
To fetch a specific version check packages page on github ([here](https://github.com/precimed/gsa-mixer/pkgs/container/gsa-mixer) for Docker, [here](https://github.com/precimed/gsa-mixer/pkgs/container/gsa-mixer_sif) for Singularity), and update the above with a specific tag, e.g. ``gsa-mixer:sha-a7b47d3``.

If you have built MiXeR's native code locally, use
```
export GSA_MIXER_ROOT=$HOME/github/precimed/gsa-mixer                   # adjust accordingly
export BGMG_SHARED_LIBRARY="$GSA_MIXER_ROOT/src/build/lib/libbgmg.so"
export MIXER_PY="python $GSA_MIXER_ROOT/precimed/mixer.py"
```

The usage of ``${MIXER_PY}`` should be the same regardless of whether you use Docker or singularity version,
however for most users we recommend running through singularity container (mainly because singularity is more commonly available in HPC clusters).
If you use docker version, you may need to customize ``$DOCKER_RUN`` variable to your environment, e.g. 
you may not need to invoke docker as sudo;
you may not need ``--user $(id -u):$(id -g)`` (this forces docker to run commands as current user);
you may also try replacing ``$PWD`` with ``pwd`` (same as in the above $MIXER_PY command using singularity container),
so that current working directory is mounted to the container even if you change it after defining ``MIXER_PY`` variable.

The above containers are only generated for CPUs with x86 architecture (e.g. intel or AMD CPUs), and do not support ARM architectures (for example the are not compatible with newer Macbook laptops with M1/M2/M3 chips).
The containers are based on the following [Dockerfile](Dockerfile), built using Github actions ([this workflow](.github/workflows/docker_build_push.yml)). We also include [scripts/from_docker_image.sh](scripts/from_docker_image.sh) shell script to convert locally built Docker container into singularity, which is only relevant if you're building these containers yourself.

### Installing Docker or Singularity on your local machine

To install Docker refer to its documentation: https://docs.docker.com/get-started/get-docker/ .
Note that by default after Docker installation you need to run it as sudo. 
The instructions below assume you can run docker as your own user.
To allow this you need to your user to 'docker' usergroup, as instructed [here](https://docs.docker.com/engine/install/linux-postinstall/).
Then check that you can run ``docker run hello-world`` command without sudo. If this works you're good to go.

Singularity software (https://sylabs.io/docs/) is most likely available in your HPC environment, however it's also
not too hard to get it up an running on your laptop (especially on Ubuntu, probably also on older MAC with an intel CPU).
To install singularity on Ubuntu follow steps described here: https://sylabs.io/guides/3.7/user-guide/quick_start.html
Note that ``sudo apt-get`` can give only a very old version of singularity, which isn't sufficient.
Therefore it's best to build singularity locally. 
Note that building singularity from source code depends on [GO](https://go.dev/doc/install), 
so it must be installed first. One you have singularity up and running, it might be usefult o have a look at
["singularity shell" options](https://sylabs.io/guides/3.2/user-guide/cli/singularity_shell.html#options) and
[Bind paths and mounts](https://sylabs.io/guides/3.2/user-guide/bind_paths_and_mounts.html) pages from the documentation.
``oras`` pre-built binary can be installed from [here](https://github.com/oras-project/oras/releases).

## Getting Started Example

This section depends on example data files located in ``precimed/mixer-test/data`` folder of this repository.
The easiest way of downloading the example data might be to just ``git clone https://github.com/precimed/gsa-mixer.git``, to clone the entire repository. Note that for this to work properly you will need ``git lfs`` to be configured (see [here](https://git-lfs.com/)).
The following files are needed to run the examples:
* ```g1000_eur_hm3_chr21to22.[bed,bim,fam]``` - EUR subset of 1kG Phase3 individuals (N=503) for M=34958 SNPs from chr21 and chr22, already constrained to HapMap3 SNPs
*  ```g1000_eur_hm3_chr[21,22].ld``` - LD matrix derived from the above genotypes using ``mixer.py ld`` command
* ```g1000_eur_hm3_chr@.snps``` - SNPs used to subset GWAS z-scores used in fit procedure; the set of SNPs is derived from the above genotypes with ``mixer.py snps`` command
* ```trait1.sumstats.gz``` and ```trait2.sumstats.gz``` - GWAS summary statistics for two traits (only the first trait is used in GSA-MiXeR demo example)
* ```partial.pheno``` - two synthesized phenotypes each with SNP-h2=0.7, generated from the above genotypes using additive genetic model; this file is not used by GSA-MiXeR, but it was used to produce the above GWAS summary statistics via ``plink2 --glm`` call.
* ```g1000_eur_hm3_chr[21,22].annot.gz``` - randomly generated functional annotations in sLDSC format
* ```go-file-baseline.csv``` - baseline model with three gene sets (all_genes, coding_genes, pseudo_genes);
* ```go-file-gene.csv``` - enrichment model with in total 435 real genes from chr21 and chr22
* ```go-file-geneset.csv``` - enrichment model with 562 real gene-sets (constrained to genes on chr21 and chr22)

GSA-MiXeR usage can be illustrated with the following steps.
The first two steps (``$MIXER ld`` and ``$MIXER snps``) are optional, as they relate to producing reference files, which for real-data analysis usually should be downloaded via the links provided below.
Note that ``@`` symbol must remain as it is in all commands, i.e. you don't need to exchange it with a specific chromosome label.
All commands below assume that demo data is locate in your current folder.
Expected execution time of all commands below on a standard laptop is less than 60 seconds.

```
cd precimed/mixer-test/data

for chri in {21..22}; do ${MIXER_PY} ld --bfile g1000_eur_hm3_chr$chri --r2min 0.05 --ldscore-r2min 0.01 --out g1000_eur_hm3_chr$chri.ld --ld-window-kb 10000; done  

${MIXER_PY} snps --bim-file g1000_eur_hm3_chr@.bim --ld-file g1000_eur_hm3_chr@.ld --chr2use 21-22 --r2 0.6 --maf 0.05 --subset 20000 --out g1000_eur_hm3_chr@.snps --seed 123

# split summary statistics into one file per chromosome
${MIXER_PY} split_sumstats --trait1-file trait1.sumstats.gz --out trait1.chr@.sumstats.gz --chr2use 21-22

# generate .bin file for --loadlib-file argument
${MIXER_PY} plsa \
      --bim-file g1000_eur_hm3_chr@.bim \
      --ld-file g1000_eur_hm3_chr@.ld \
      --use-complete-tag-indices --chr2use 21-22 --exclude-ranges [] \
      --savelib-file g1000_eur_hm3_chr@.bin \
      --out g1000_eur_hm3_chr@

# fit baseline model, and use it to calculate heritability attributed to gene-sets in go-file-geneset.csv
${MIXER_PY} plsa --gsa-base \
--trait1-file trait1.chr@.sumstats.gz \
--use-complete-tag-indices \
--bim-file g1000_eur_hm3_chr@.bim \
--loadlib-file g1000_eur_hm3_chr@.bin \
--annot-file g1000_eur_hm3_chr@.annot.gz \
--go-file go-file-baseline.csv \
--extract g1000_eur_hm3_chr@.snps \
--exclude-ranges chr21:20-21MB chr22:19100-19900KB \
--chr2use 21-22 --seed 123 \
--adam-epoch 3 3 --adam-step 0.064 0.032 \
--out plsa_base

# fit enrichment model, and use it to calculate heritability attributed to gene-sets in go-file-geneset.csv
${MIXER_PY} plsa --gsa-full \
--trait1-file trait1.chr@.sumstats.gz \
--use-complete-tag-indices \
--bim-file g1000_eur_hm3_chr@.bim \
--loadlib-file g1000_eur_hm3_chr@.bin \
--annot-file g1000_eur_hm3_chr@.annot.gz \
--go-file go-file-gene.csv \
--go-file-test go-file-geneset.csv \
--extract g1000_eur_hm3_chr@.snps \
--load-params-file plsa_base.json \
--exclude-ranges chr21:20-21MB chr22:19100-19900KB \
--chr2use 21-22 --seed 123 \
--adam-epoch 3 3 --adam-step 0.064 0.032 \
--out plsa_full
```

The commands above are customized to run the analysis faster.
For real-data analysis the commands will need to be adjusted. 
The [scripts/GSA_MIXER.job](scripts/GSA_MIXER.job) is a good starting point for real-world example of the GSA-MiXeR application; note how this script implements the following changes, as compared to the above commands from the getting started example:
* remove ``--exclude-ranges chr21:20-21MB chr22:19100-19900KB``; by default ``--exclude-ranges`` will exclude MHC region
* remove ``--chr2use 21-22``; by default ``--chr2use`` applies to all chromosomes
* remove ``--adam-epoch 3 3 --adam-step 0.064 0.032``, as this stops Adam fit procedure too early
* remove ``--extract g1000_eur_hm3_chr@.snps``, to use all available SNPs for fitting the model.
The multi-start procedure involving 20 re-runs of the fit procedure, with each constrainted to a random subset of SNPs, was only relevant to cross-trait MiXeR and is not recommended for GSA-MiXeR.

## Input Data Formats

GSA-MiXeR format for summary statistics (``--trait1-file``) is compatible with LD Score Regression
(i.e. the ``.sumstats.gz`` files), and must include the following columns:
* Eithe one of the following:
  * ``SNP`` or ``RSID`` (marker name), or
  * ``CHR`` (chromosome label) and  ``BP`` or ``POS`` (genomic corrdinates), in a build that is compatible with the reference build (``--bim-file`` argument)
* ``A1`` or ``EffectAllele`` (reference allele)
* ``A2`` or ``OtherAllele`` (alternative allele)
* ``N`` (sample size); for case-control studies this should be the effective sample size computed as ``N=4/(1/ncases+1/ncontrols)``
* ``Z`` (signed test statistic)
Column names must be exactly as defined above, except for upper/lower case which can be arbitrary (all column names from the input file are converted to lower case prior to matching them with expected column names defined above).

It's beneficial to have both ``SNP`` and ``CHR``/``BP`` columns in the data.
In this situation matching SNPs with the reference (``--bim-file``) is performed on marker name (``SNP`` column).
For the remaining SNPs GSA-MiXeR attempts to match using CHR:BP:A1:A2 codes, accounting for situations when (1) ``A1``/``A2`` are swapped between GWAS summary statistics and the reference, (2) alleles are coded on different strands, (3) both of the above. 
The sign of z-score is refersed when needed during this procedure.

We advice against filtering down summary statistics to the set of SNPs included in HapMap3.
Prior to running GSA-MiXeR we advice filtering SNPs with bad imputation quality, if ``INFO`` column is available in summary statistics.
If per-SNP sample size is available, we advice filtering out SNPs with N below half of the median sample size across SNPs.
Other filtering options are built into GSA-MiXeR software, including ``--exclude-ranges`` option to filter out special regions such as MHC, and ``--maf`` and ``--randprune-maf`` to filter out based on minor allele frequency.

Prior to running GSA-MiXeR you will need to split summary statistics into one file per chromosome, as shown in [GSA_MIXER.job](scripts/GSA_MIXER.job):
```
${MIXER_PY} split_sumstats \
    --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats.gz
    --out ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.chr@.sumstats.gz
```

## Download LD matrix and other reference files

All reference data described below is based on EUR ancestry, and use ``hg19`` / ``GRCh37`` genomic build.

Reference files derived from 1kG Phase3 EUR population are available for download from [here](https://github.com/comorment/mixer/tree/main/reference/ldsc/1000G_EUR_Phase3_plink).
The easiest way of downloading the example data might be to just ``git clone https://github.com/comorment/mixer.git``, to clone the entire [comorment/mixer](https://github.com/comorment/mixer) repository. Note that for this to work properly you will need ``git lfs`` to be configured (see [here](https://git-lfs.com/)).
The following files are needed:

```
1000G_EUR_Phase3_plink/1000G.EUR.QC.[1-22].bim                      # ``--bim-file`` argument
1000G_EUR_Phase3_plink/baseline_v2.2_1000G.EUR.QC.[1-22].annot.gz   # ``--annot-file`` / ``--annot-file-test`` arguments
1000G_EUR_Phase3_plink/1000G.EUR.QC.[1-22].run4.ld                  # ``--ld-file`` argument
1000G_EUR_Phase3_plink/1000G.EUR.QC.@.[bed/bim,fam]                 # reference for MAGMA analysis, merged across chromosomes
```

Functional annotations are derived from [sLDSC baselineLD_v2.2](https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/baselineLD_v2.2_bedfiles.tgz) using scripts from [here](https://github.com/ofrei/eas_partitioned_ldscore) to annotate UKB, HRC and 1kG  references. Note that one does not need to compute LD-scores for these annotations, because MiXeR does this internally using sparse LD matrix stored in ``--ld-file`` it receives as an argument.


Additionally you need to download gene- and gene-set definitions from [here](https://github.com/precimed/gsa-mixer/tree/main/reference/) (for future reference, these definitions are derived using scripts from [here](https://github.com/ofrei/genesets/blob/main/prepare_genesets_v2.ipynb) ):
```
gsa-mixer-baseline-annot_10mar2023.csv           # ``--go-file`` (baseline model)
gsa-mixer-gene-annot_10mar2023.csv               # ``--go-file`` (model)
gsa-mixer-geneset-annot_10mar2023.csv            # ``--go-file-test`` (only gene-sets)
gsa-mixer-genesetLOO-annot_10mar2023.csv         # ``--go-file-test`  (only gene-sets, with leave-one-gene-out)
gsa-mixer-hybrid-annot_10mar2023.csv             # ``--go-file-test`  (genes and gene-sets)
gsa-mixer-hybridLOO-annot_10mar2023.csv          # ``--go-file-test`  (genes and gene-sets, with leave-one-gene-out)
magma-gene-annot_10mar2023.csv                   # gsa-mixer-gene-annot_10mar2023.csv converted to MAGMA format
magma-geneset-annot_10mar2023.csv                # gsa-mixer-geneset-annot_10mar2023.csv converted to MAGMA format
```

After downloading LD reference files we advice using ``--savelib-file`` option as shown below to generate ``.bin`` files,
with compressed representation of the reference. After that loading reference is possible with ``--loadlib-file``, providing considerable speedup over passing ``--ld-file`` argument.

The reference needs to be saved in two formats. The following example produces ``.bin`` files for ``plsa`` analysis, yielding its own ``.bin`` file for each chromosome. The ``--savelib-file`` argument must include ``@`` symbol which will be replaced with an actual chromosome label.
In order to reduce peak memory use the command is executed through a ``for`` loop, i.e. separately for each chromosome, but it's also ok to remove ``--chr2use $chr`` option and run just a single `${MIXER_PY} plsa --savelib-file` command.
```
for chri in {1..22}; do ${MIXER_PY} plsa \
      --bim-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --use-complete-tag-indices --exclude-ranges [] --chr2use $chri \
      --savelib-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin \
      --out 1000G_EUR_Phase3_plink/1000G.EUR.QC.@; done
```

The following example produces ``.bin`` file for ``fit1``,``fit2``,``test1``,``test2`` steps (cross-trait MiXeR; not relevant for GSA-MiXeR), yielding a single ``.bin`` file combined across all chromosomes. This has a fairly high peak memory usage. The ``--savelib-file`` argument does not need to include ``@`` symbol (if it does, the ``@`` symbol will stay unchanged, and simply be part of the output file name):
```
${MIXER_PY} fit1 \
      --bim-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --use-complete-tag-indices --exclude-ranges [] \
      --savelib-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin \
      --out 1000G_EUR_Phase3_plink/1000G.EUR.QC.@
```

## Perform GSA-MiXeR and MAGMA analyses

[scripts/GSA_MIXER.job](scripts/GSA_MIXER.job) script is a good starting point for performing GSA-MiXeR and MAGMA analyses. Below is an overview of its key sections.

### Define SLURM parameters

GSA-MiXeR analysis takes around 6 to 12 hours using 8-core machine, and utilize around 30 GB of RAM. The following configuration might be a good starting point:
```
#SBATCH --job-name=gsamixer
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account=nn9114k
#SBATCH --mem-per-cpu=4569M   # 178.5 GB / 40 cores - https://documentation.sigma2.no/jobs/job_types/saga_job_types.html
#SBATCH --cpus-per-task=8     # if you change keep --threads $THREADS argument in sync
```

### Define data location, and singularity-related variables
```
export GITHUB=/cluster/projects/nn9114k/github
export MIXER_SIF=${GITHUB}/precimed/gsa-mixer/containers/latest/gsa-mixer.sif
export SUMSTATS_FOLDER=/cluster/projects/nn9114k/oleksanf/gsa-mixer/sumstats
export SUMSTATS_FILE=PGC_SCZ_0518_EUR
export OUT_FOLDER=/cluster/projects/nn9114k/oleksanf/gsa-mixer/out2
export BIND="--bind /cluster/projects/nn9114k:/cluster/projects/nn9114k"

export REFERENCE_FOLDER=${GITHUB}/precimed/gsa-mixer/reference
export BIM_FILE=${GITHUB}/comorment/mixer/reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim
export LOADLIB_FILE=${GITHUB}/comorment/mixer/reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin
export ANNOT_FILE=${GITHUB}/comorment/mixer/reference/ldsc/1000G_EUR_Phase3_plink/baseline_v2.2_1000G.EUR.QC.@.annot.gz

export MAGMA_BFILE=${GITHUB}/comorment/mixer/reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@
export MAGMA_GENE_LOC=${GITHUB}/precimed/gsa-mixer/reference/magma-gene-annot_10mar2023.csv
export MAGMA_SET_ANNOT=${GITHUB}/precimed/gsa-mixer/reference/magma-geneset-annot_10mar2023.csv

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

### Perform MAGMA analysis:

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

## Command-line reference

```
usage: mixer.py plsa [-h] [--out OUT] [--lib LIB] [--log LOG]
                     [--bim-file BIM_FILE] [--ld-file LD_FILE]
                     [--chr2use CHR2USE] [--extract EXTRACT]
                     [--exclude EXCLUDE]
                     [--exclude-ranges EXCLUDE_RANGES [EXCLUDE_RANGES ...]]
                     [--allow-ambiguous-snps] [--trait1-file TRAIT1_FILE]
                     [--z1max Z1MAX] [--randprune-n RANDPRUNE_N]
                     [--randprune-r2 RANDPRUNE_R2]
                     [--randprune-maf RANDPRUNE_MAF] [--seed SEED]
                     [--threads THREADS [THREADS ...]]
                     [--kmax KMAX [KMAX ...]] [--annot-file ANNOT_FILE]
                     [--annot-file-test ANNOT_FILE_TEST] [--go-file GO_FILE]
                     [--go-file-test GO_FILE_TEST]
                     [--go-all-genes-label GO_ALL_GENES_LABEL]
                     [--go-extend-bp GO_EXTEND_BP] [--sig2-zeroL SIG2_ZEROL]
                     [--s-value S_VALUE] [--l-value L_VALUE]
                     [--pi-value PI_VALUE] [--s-init S_INIT] [--l-init L_INIT]
                     [--pi-init PI_INIT] [--annot-p ANNOT_P] [--gene-p GENE_P]
                     [--adam-epoch ADAM_EPOCH [ADAM_EPOCH ...]]
                     [--adam-beta1 ADAM_BETA1] [--adam-beta2 ADAM_BETA2]
                     [--adam-eps ADAM_EPS]
                     [--adam-step ADAM_STEP [ADAM_STEP ...]] [--adam-disable]
                     [--load-params-file LOAD_PARAMS_FILE] [--make-snps-file]

optional arguments:
  -h, --help            show this help message and exit
  --out OUT             prefix for the output files, such as <out>.json
                        (default: mixer);
  --lib LIB             path to libbgmg.so plugin (default: libbgmg.so); can
                        be also specified via BGMG_SHARED_LIBRARY env
                        variable.
  --log LOG             file to output log (default: <out>.log); NB! if --log
                        points to an existing file the new lines will be
                        appended to it at the end of the file.
  --bim-file BIM_FILE   plink bim file (required argument); defines the
                        reference set of SNPs used for the analysis. Marker
                        names must not have duplicated entries. May contain
                        symbol '@', which will be replaced by an actual
                        chromosome label.
  --ld-file LD_FILE     file with linkage disequilibrium information,
                        generated via 'mixer.py ld' command (required
                        argument); may contain symbol '@', similarly to --bim-
                        file argument.
  --chr2use CHR2USE     chromosome labels to use (default: 1-22); chromosome
                        must be labeled by integer, i.e. X and Y are not
                        acceptable; example of valid arguments: '1,2,3' or
                        '1-4,12,16-20'
  --extract EXTRACT     (optional) File with variants to include in the
                        analysis. By default, all variants are included. This
                        applies to GWAS tag SNPs, however LD is still computed
                        towards the entire reference provided by --bim-file.
                        This applies before --exclude, so a variant listed in
                        --exclude won't be re-introduced if it's also present
                        in --extract list.
  --exclude EXCLUDE     (optional) File with variants to exclude from the
                        analysis.
  --exclude-ranges EXCLUDE_RANGES [EXCLUDE_RANGES ...]
                        (default: ['MHC']) exclude SNPs in ranges of base pair
                        position; the syntax is chr:from-to, for example
                        6:25000000-35000000; multiple regions can be excluded;
                        "chr" prefix prior to chromosome label, as well as KB
                        and MB suffices are allowed, e.g. chr6:25-35MB is a
                        valid exclusion range. Some special case regions are
                        also supported, for example "--exclude-ranges MHC
                        APOE".To overwrite the defaul, pass "--exclude ranges
                        []".
  --allow-ambiguous-snps
                        advanced option (expert use only); a flag allowing to
                        include A/T and C/G SNPs in fit procedure.
  --trait1-file TRAIT1_FILE
                        GWAS summary statistics for the first trait (required
                        argument); for 'plsa' analysis it is recommended to
                        split GWAS summary statistics per chromosome; if this
                        is done then --trait1-file should contain symbol '@',
                        which will be replaced by an actual chromosome label.
  --z1max Z1MAX         right-censoring threshold for the first trait
                        (default: None); recommended setting: '--z1max 9.336'
                        (equivalent to p-value 1e-20)
  --randprune-n RANDPRUNE_N
                        number of random pruning iterations (default: 64)
  --randprune-r2 RANDPRUNE_R2
                        threshold for random pruning (default: 0.1)
  --randprune-maf RANDPRUNE_MAF
                        threshold for minor allele frequency (default: 0.05);
                        applies to tag SNPs to include in the fit procedure
  --seed SEED           random seed (default: None).
  --threads THREADS [THREADS ...]
                        specify how many concurrent threads to use for
                        computations; (default: total number of CPU cores
                        available)
  --kmax KMAX [KMAX ...]
                        number of sampling iterations for log-likelihod and
                        posterior delta (default: 20000)
  --annot-file ANNOT_FILE
                        (optional) path to binary annotations in LD score
                        regression format, i.e. <path>/baseline.@.annot.gz for
                        fitting enrichment model model. This must include the
                        first column with all ones ('base' annotation category
                        covering the entire genome). Enrichment scores
                        computed for --annot-file will be saved to
                        <out>.annot_enrich.csv
  --annot-file-test ANNOT_FILE_TEST
                        (optional) path to binary annotations in LD score
                        regression format, i.e. <path>/baseline.@.annot.gz for
                        evaluating enrichment. If provided, enrichment scores
                        computed for --annot-file-test will be saved to
                        <out>.annot_test_enrich.csv.
  --go-file GO_FILE     (optional) path to GO antology file for fitting
                        enrichment model model. The format is described in the
                        documentation. 'base' category that covers entire
                        genome will be added automatically. Enrichment scores
                        computed for --go-file will be saved to
                        <out>.go_enrich.csv
  --go-file-test GO_FILE_TEST
                        (optional) path to an additional GO antology file for
                        evaluating enrichment, same convention as --go-file.
                        If provided, enrichment scores computed for --go-test-
                        file will be saved to <out>.go_test_enrich.csv
  --go-all-genes-label GO_ALL_GENES_LABEL
                        reference gene-set to calibrate fold enrichment, e.g.
                        allowing to compute enrichment w.r.t. the set of all
                        coding genes (default:'base')
  --go-extend-bp GO_EXTEND_BP
                        extends each gene by this many base pairs, defining a
                        symmetric window up and downstream (default: 10000)
  --sig2-zeroL SIG2_ZEROL
                        (optional) constraint for 'sig2_zeroL' parameter;
                        recommended setting for gene-set enrichment analysis:
                        '--sig2-zeroL 0'
  --s-value S_VALUE     (optional) constraint for the 's' parameter;
                        recommended setting for gene-set enrichment analysis:
                        '--s-value -0.25'
  --l-value L_VALUE     (optional) constraint for the 'l' parameter
  --pi-value PI_VALUE   (optional) constraint for the 'pi' parameter;
                        recommended setting for gene-set enrichment analysis:
                        '--pi-value 1.0'
  --s-init S_INIT       initial value for the 's' parameter (default: -0.25);
                        does not apply if --s-value is provided)
  --l-init L_INIT       initial value for the 'l' parameter (default: -0.125);
                        does not apply if --l-value is provided)
  --pi-init PI_INIT     initial value for the 'pi' parameter (default: 0.001);
                        does not apply if --pi-value is provided)
  --annot-p ANNOT_P     power factor for sigma2 aggregation in overlapping
                        annotations (default: 1)
  --gene-p GENE_P       power factor for sigma2 aggregation in overlapping
                        gene-sets (default: 1)
  --adam-epoch ADAM_EPOCH [ADAM_EPOCH ...]
                        number of iterations in ADAM procedure (default: [10,
                        10, 10, 10, 10, 10, 10, 10, 10, 10])
  --adam-beta1 ADAM_BETA1
                        beta_1 parameter in ADAM procedure (default: 0.9)
  --adam-beta2 ADAM_BETA2
                        beta_2 parameter in ADAM procedure (default: 0.99)
  --adam-eps ADAM_EPS   epsilon parameter in ADAM procedure (default: 1e-08)
  --adam-step ADAM_STEP [ADAM_STEP ...]
                        step parameter in ADAM procedure (default: [0.064,
                        0.032, 0.016, 0.008, 0.004, 0.002, 0.001, 0.0005,
                        0.00025, 0.0001])
  --adam-disable        a flag allowing to disable optimization; typical
                        usecase would be in conjunction with these flags: '--
                        adam-disable --load-params-file <out-of-a-previous-
                        run>.json --make-snps-file --allow-ambiguous-snps'
  --load-params-file LOAD_PARAMS_FILE
                        (optional) params of the fitted model to load;
                        expected to be from a 'mixer.py plsa' run
  --make-snps-file      a flag allowing to generate file with per-SNP
                        estimates; will generate <out>.snps.csv output file
```

## Generate LD Reference

Analyses in the original publication are based on UKB and HRC reference panen, partly shared here:
* [UKB reference](https://github.com/precimed/gas-mixer/tree/main/reference/ukb_EUR_qc) and 
* [HRC reference](https://github.com/precimed/gsa-mixer/tree/main/reference/hrc_EUR_qc).

These references require the following files:
```
ukb_EUR_qc/about_UKB_qc.txt                                 # overview of QC procedure
ukb_EUR_qc/ukb_imp_chr[1-22]_v3_qc.bim                      # ``--bim-file`` argument
ukb_EUR_qc/baseline_v2.2_ukb_imp_chr[1-22]_v3_qc.annot.gz   # ``--annot-file`` / ``--annot-file-test`` arguments
ukb_EUR_qc/ukb_imp_chr[1-22]_v3_qc.run1.ld                  # ``--ld-file`` argument (files not shared)
ukb_EUR_qc/ukb_imp_chr@_qc.prune_rand2M_rep[1-20].snps      # ``--extract`` argument

hrc_EUR_qc/about_HRC_qc.txt                                 # overview of QC procedure
hrc_EUR_qc/hrc_chr[1-22]_EUR_qc.bim                         # ``--bim-file`` argument
hrc_EUR_qc/hrc_chr[1-22]_EUR_qc.run1.ld                     # ``--annot-file`` / ``--annot-file-test`` arguments
hrc_EUR_qc/baseline_v2.2_hrc_chr[1-22]_EUR_qc.annot.gz      # ``--ld-file`` argument (files not shared)
hrc_EUR_qc/hrc_EUR_qc.prune_rand2M_rep[1-20].snps           # ``--extract`` argument
```
All files are shared except for the full LD matrix, which is not shared due to concerns of de-identifying the subjects.
To re-generate LD matrices you will need to obtain access to individual-level data of UKB or HRC subjects, and re-run analysis as defined in ``about_UKB_qc.txt`` / ``about_HRC_qc.txt`` steps.

GSA-MiXeR reference files can be prepared from plink bfile using ``mixer.py ld`` and ``mixer.py snps`` commands as described below. It's important that the reference genotypes contain unrelated individuals only, constrained to a single population.
Note that ``@`` symbol must remain as it is in all commands, i.e. you don't need to exchange it with a specific chromosome label.

Compute LD matrices (one per chromosome), later to be used with ``--ld-file`` argument.
```
#!/bin/bash
#SBATCH --job-name=gsamixer
#SBATCH --account=p697
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22

export CHR=${SLURM_ARRAY_TASK_ID}
export MIXER_SIF=mixer.sif
export MIXER_PY="singularity exec --home pwd:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

${MIXER_PY} ld --bfile chr${CHR} --r2min 0.01 --ldscore-r2min 0.0001 --ld-window-kb 10000 --out chr${CHR}.ld
```

Compute SNPs subsets (one for each of 20 random repeats), later to be used with ``--extract`` argument.
This step is still relevant for cross-trait MiXeR, but this is not used by GSA-MiXeR.
If you generate custom reference for GSA-MiXeR, this step can be skiped.
```
#!/bin/bash
#SBATCH --job-name=gsamixer
#SBATCH --account=p697
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=8
#SBATCH --array=1-20

export REP=${SLURM_ARRAY_TASK_ID}
export MIXER_SIF=gsa-mixer.sif
export MIXER_PY="singularity exec --home pwd:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

${MIXER_PY} snps --bim-file chr${CHR} --ld-file chr@ --chr2use 1-22 --maf 0.05 --subset 3000000 --seed $REP --out rep${REP}.snps
```

If only a random subset of SNPs is needed it's faster to use linux's ``cut`` and ``shuf`` commands:
```
for i in {1..20}
do 
  cat hrc_chr*_EUR_qc.bim | cut -f 2 | shuf | head -n 2000000 | sort > hrc_EUR_qc.prune_rand2M_rep$i.snps
done
```

Full command-line reference for ``mixer.py ld`` and ``mixer.py snps`` is as follows:
```
usage: mixer.py ld [-h] [--out OUT] [--lib LIB] [--log LOG] [--bfile BFILE]
                   [--r2min R2MIN] [--ldscore-r2min LDSCORE_R2MIN]
                   [--ld-window-kb LD_WINDOW_KB] [--ld-window LD_WINDOW]

optional arguments:
  -h, --help            show this help message and exit
  --out OUT             prefix for the output files, such as <out>.json
                        (default: mixer);
  --lib LIB             path to libbgmg.so plugin (default: libbgmg.so); can
                        be also specified via BGMG_SHARED_LIBRARY env
                        variable.
  --log LOG             file to output log (default: <out>.log); NB! if --log
                        points to an existing file the new lines will be
                        appended to it at the end of the file.
  --bfile BFILE         path to plink bfile (required argument)
  --r2min R2MIN         r2 values above this threshold will be stored in
                        sparse LD format (default: 0.05)
  --ldscore-r2min LDSCORE_R2MIN
                        r2 values above this threshold (and below --r2min)
                        will be stored as LD scores that contribute to the
                        cost function via an infinitesimal model (default:
                        0.001)
  --ld-window-kb LD_WINDOW_KB
                        limit window similar to --ld-window-kb in 'plink r2';
                        0 will disable this constraint (default: 0); either
                        ld-window-kb or ld-window argument must be provided
  --ld-window LD_WINDOW
                        limit window similar to --ld-window in 'plink r2'; 0
                        will disable this constraint (default: 0); either ld-
                        window-kb or ld-window argument must be provided
```

```
usage: mixer.py snps [-h] [--out OUT] [--lib LIB] [--log LOG]
                     [--bim-file BIM_FILE] [--ld-file LD_FILE]
                     [--chr2use CHR2USE] [--r2 R2] [--maf MAF]
                     [--subset SUBSET] [--seed SEED]

optional arguments:
  -h, --help           show this help message and exit
  --out OUT            prefix for the output files, such as <out>.json
                       (default: mixer);
  --lib LIB            path to libbgmg.so plugin (default: libbgmg.so); can be
                       also specified via BGMG_SHARED_LIBRARY env variable.
  --log LOG            file to output log (default: <out>.log); NB! if --log
                       points to an existing file the new lines will be
                       appended to it at the end of the file.
  --bim-file BIM_FILE  plink bim file (required argument); defines the
                       reference set of SNPs used for the analysis. Marker
                       names must not have duplicated entries. May contain
                       symbol '@', which will be replaced by an actual
                       chromosome label.
  --ld-file LD_FILE    file with linkage disequilibrium information, generated
                       via 'mixer.py ld' command (required argument); may
                       contain symbol '@', similarly to --bim-file argument.
  --chr2use CHR2USE    chromosome labels to use (default: 1-22); chromosome
                       must be labeled by integer, i.e. X and Y are not
                       acceptable; example of valid arguments: '1,2,3' or
                       '1-4,12,16-20'
  --r2 R2              LD r2 threshold for prunning SNPs (default: 0.8)
  --maf MAF            minor allele frequence (MAF) threshold (default: 0.05)
  --subset SUBSET      number of SNPs to randomly select (default: 2000000)
  --seed SEED          random seed (default: None)
```
