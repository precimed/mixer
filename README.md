## User documentation for MiXeR analyses (univariate, bivariate and GSA-MiXeR)

This repository (https://github.com/precimed/mixer) provides user documentation for MiXeR analyses (univariate, bivariate and GSA-MiXeR).
You should use this repository for most up to date instructions on how to install and run MiXeR.
Another relevant repository is https://github.com/comorment/mixer where you will find reference data.
These two repositories contain everything you need to run MiXeR.

If you previously ran univariate or cross-trait MiXeR using ``mixer.sif`` container obtained from https://github.com/comorment/mixer project, the only step you need to do to upgrade to the latest version is to re-download ``mixer.sif`` container as described [here](#install). Everything else (input data, scripts performing the analysis, and output formats) remain the same, and the results should be nearly identical, except for minor bug fixes.

If you ancounter an issue, please submit a ticket: https://github.com/precimed/mixer/issues/new . In the past my response to new tickets was very bad. If you've already submitted a ticket but didn't get a response please give me a second chance to address your issue - just push a comment and tag me @ofrei , if your question is still relevant .

MiXeR source code is in https://github.com/precimed/gsa-mixer repository. This also covers univariate and cross-trait analyses, i.e. not just GSA MiXeR.
Either way, it is only relevant to developers interested to contribute pull requests to mixer code.

Additional instructions for users at the NORMENT centre are available in https://github.com/precimed/mixer_at_norment.
If the link doesn't work please reach me out by e-mail to get access.

Kind regards,
Oleksandr Frei.

## Contents

* [Introduction](#introduction)
* [Install procedures](#install)
* [Getting started](#getting-started)
  * [GSA MiXeR hello-world example](#gsa-mixer-hello-world-example)
  * [Cross-trait MiXeR hello-world example](#cross-trait-mixer-hello-world-example)
* [Tutorials (examples of real dat analyses)](#tutorials)
* [GWAS summary statistics format](#gwas-summary-statistics-format)
* [Reference data - downloads](#reference-data---downloads)
  * [Generate your own LD reference](#generate-your-own-ld-reference-cross-trait--gsa)
* [FAQ](#faq)
  * [AIC BIC interpretation](#aic-bic-interpretation) 
  * [Upgrade nodes from MiXeR v1.2 to v1.3](#upgrade-nodes-from-mixer-v12-to-v13)
  * [Installing Docker, Singularity, oras on your local machine](#installing-docker-singularity-oras-on-your-local-machine)
* [Changelog](#changelog)

## Introduction

**GSA-MiXeR** is a new technique for competitive gene-set analysis, which fits a model for gene-set heritability enrichments for complex human traits, thus allowing the quantification of partitioned heritability and fold enrichment for small gene-sets.

**Cross-trait MiXeR** is a statistical tool which quantifies polygenic overlap between complex traits irrespective of genetic correlation, using GWAS summary statistics. MiXeR results are presented as a Venn diagram of unique and shared polygenic components across traits. Cross-trait (bivariate) analysis builds on univariate analysis quantifying polygenicity of the traits.

If you use this package, please cite the original work and all reference data used in your analysis.

* for univariate analysis: D. Holland et al., Beyond SNP Heritability: Polygenicity and Discoverability Estimated for Multiple Phenotypes with a Univariate Gaussian Mixture Model, PLOS Genetics, 2020, https://doi.org/10.1371/journal.pgen.1008612
* for cross-trait analysis: O.Frei et al., Bivariate causal mixture model quantifies polygenic overlap between complex traits beyond genetic correlation, Nature Communications, 2019, https://www.nature.com/articles/s41467-019-10310-0
* for gene-set enrichment analysis: O.Frei et al., Improved functional mapping with GSA-MiXeR implicates biologically specific gene-sets and estimates enrichment magnitude, Nature Genetics, 2024, https://www.nature.com/articles/s41588-024-01771-1

## Install

MiXeR software is released both as a Docker container, and as a pre-compiled singularity (apptainer) container. Use the following commands to check if you have Docker and/or singulrity available in your environment:
```
# check if Docker software is installed
>docker --version
Docker version 20.10.7, build 20.10.7-0ubuntu5~21.04.2

# check if singularity software is installed
>singularity --version
singularity version 3.7.4
```

Check [this page](https://github.com/precimed/gsa-mixer/releases) to find tag for the latest release, e.g. ``v2.2.1``.
Don't use versions marked as pre-release.

Then, to dowload Docker version of the GSA-MiXeR, use the following command:
```
docker pull ghcr.io/precimed/gsa-mixer:2.2.1   # replace 2.2.1 with latest tag without "v"
export DOCKER_RUN="docker run -v $PWD:/home -w /home"
export MIXER_PY="$DOCKER_RUN ghcr.io/precimed/gsa-mixer:latest python /tools/mixer/precimed/mixer.py"
```
To download singularity version of the GSA-MiXeR, use the following command:
```
oras pull ghcr.io/precimed/gsa-mixer_sif:2.2.1 && mv gsa-mixer.sif mixer.sif
export MIXER_SIF=<path>/mixer.sif
export MIXER_PY="singularity exec --home pwd:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"
```

To test your installation, try
```
${MIXER_PY} --version
${MIXER_PY} --help
```

The usage of ``${MIXER_PY}`` should be the same regardless of whether you use Docker or singularity version,
however for most users we recommend running through singularity container (mainly because singularity is more commonly available in HPC clusters).
If you use docker version, you may need to customize ``$DOCKER_RUN`` variable to your environment, e.g. 
you may also try replacing ``$PWD`` with ``pwd`` so that current working directory is correctly mounted to the container even if you change it after defining ``MIXER_PY`` variable.

> [!WARNING] 
> The above containers are only generated for CPUs with x86 architecture (e.g. intel or AMD CPUs), and do not support ARM architectures (for example the are not compatible with newer Macbook laptops with M1/M2/M3 chips).

## Getting started

Toy examples below depend only on data files from [this folder](https://github.com/comorment/mixer/tree/main/reference/mixer_hello_world). You can download them as follows without downloading full reference files:
```
wget https://raw.githubusercontent.com/comorment/mixer/refs/heads/main/reference/mixer_hello_world/mixer_hello_world.tar.gz
tar -xzvf mixer_hello_world.tar.gz
```

### GSA MiXeR (hello-world example)

The following files are needed to run the GSA-MiXeR hello-world example:
* ```g1000_eur_hm3_chr@.[bed,bim,fam]``` - EUR subset of 1kG Phase3 individuals (N=503) for M=34958 SNPs from chr21 and chr22, already constrained to HapMap3 SNPs
* ```trait1.sumstats.gz``` and ```trait2.sumstats.gz``` - GWAS summary statistics for two traits (only the first trait is used in GSA-MiXeR demo example)
* ```g1000_eur_hm3_chr[21,22].annot.gz``` - randomly generated functional annotations in sLDSC format
* ```go-file-baseline.csv``` - baseline model with three gene sets (all_genes, coding_genes, pseudo_genes);
* ```go-file-gene.csv``` - enrichment model with in total 435 real genes from chr21 and chr22
* ```go-file-geneset.csv``` - enrichment model with 562 real gene-sets (constrained to genes on chr21 and chr22)

The first two steps (``$MIXER ld`` and ``$MIXER snps``) are optional, as they relate to producing reference files, which for real-data pre-generated and available for download.
Note that ``@`` symbol must remain as it is in all commands, i.e. you don't need to exchange it with a specific chromosome label.
All commands below assume that demo data is locate in your current folder.
Expected execution time of all commands below on a standard laptop is less than 60 seconds.

```
for chri in {21..22}; do ${MIXER_PY} ld --bfile g1000_eur_hm3_chr$chri --r2min 0.05 --ldscore-r2min 0.01 --out g1000_eur_hm3_chr$chri.ld --ld-window-kb 10000; done  

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
The multi-start procedure involving 20 re-runs of the fit procedure, with each constrainted to a random subset of SNPs, was only relevant to cross-trait MiXeR and is not recommended for GSA-MiXeR.

### Cross-trait MiXeR (hello-world example)

Using the same setup as above one can run the following commands to try out cross-trait analysis on a dummy data:

```
${MIXER_PY} snps --bim-file g1000_eur_hm3_chr@.bim --ld-file g1000_eur_hm3_chr@.ld --chr2use 21-22 --r2 0.6 --maf 0.05 --subset 20000 --out g1000_eur_hm3_chr@.snps --seed 123

export MIXER_FASTRUN_ARGS="--fit-sequence diffevo-fast neldermead-fast --diffevo-fast-repeats 2 "
export MIXER_COMMON_ARGS="--chr2use 21-22 --seed 123  --ld-file g1000_eur_hm3_chr@.ld --bim-file g1000_eur_hm3_chr@.bim --extract g1000_eur_hm3_chr@.snps --kmax-pdf 10 --downsample-factor 1000 "

${MIXER_PY} fit1 $MIXER_COMMON_ARGS $MIXER_FASTRUN_ARGS --trait1-file trait1.sumstats.gz --out trait1.fit
${MIXER_PY} fit1 $MIXER_COMMON_ARGS $MIXER_FASTRUN_ARGS --trait1-file trait2.sumstats.gz --out trait2.fit
${MIXER_PY} fit2 $MIXER_COMMON_ARGS $MIXER_FASTRUN_ARGS --trait1-file trait1.sumstats.gz --trait2-file trait2.sumstats.gz --trait1-params trait1.fit.json --trait2-params trait2.fit.json --out trait1_vs_trait2.fit

${MIXER_PY} test1 $MIXER_COMMON_ARGS --trait1-file trait1.sumstats.gz --load-params trait1.fit.json --out trait1.test
${MIXER_PY} test1 $MIXER_COMMON_ARGS --trait1-file trait2.sumstats.gz --load-params trait2.fit.json --out trait2.test
${MIXER_PY} test2 $MIXER_COMMON_ARGS --trait1-file trait1.sumstats.gz --trait2-file trait2.sumstats.gz --load-params trait1_vs_trait2.fit.json --out trait1_vs_trait2.test
```

The commands above are customized to run the analysis faster.
For real-data analysis the commands will need to be adjusted. 
The [usecases/mixer_real/MIXER_REAL.job](usecases/mixer_real/MIXER_REAL.job) is a good starting point for real-world example of the cross-trait MiXeR application; note how this script implements the following changes, as compared to the above commands from the getting started example:
* remove ``--chr2use 21-22``; by default ``--chr2use`` applies to all chromosomes
* remove ``--fit-sequence diffevo-fast neldermead-fast --diffevo-fast-repeats 2`` flags
* remove ``--kmax-pdf 10 --downsample-factor 1000`` flags

## Tutorials

These tutorials depend on pre-downloaded reference data, as described [here](#reference-data---downloads)

* See [usecases/mixer_simu.md](usecases/mixer_simu.md) for cross-trait MiXeR analysis on a syntetic data, 
show-casing unique trait architecture vs partial and full polygenic overlap scenarios.

* See [usecases/mixer_real.md](usecases/mixer_real.md) for tutorial on cross-trait MiXeR analysis
More information is in [the legacy README file](CROSS_TRAIT_README.md) describing cross-trait analysis.

* See jupyter notebook [usecases/run_mixer.ipynb](usecases/run_mixer.ipynb) for submitting a bunch of MiXeR jobs at once, with lists of primary and secondary traits.

* See [usecases/cross_trait_legacy_tutorial.md](usecases/cross_trait_legacy_tutorial.md) for legacy tutorial using cross-trait MiXeR.

* See [usecases/gsa_mixer.md](usecases/gsa_mixer.md) for an overview of GSA-MiXeR analysis.

## GWAS summary statistics format

Input data for MiXeR consists of summary statistics from a GWAS, and a reference panel.
Format for summary statistics (``--trait1-file``) is compatible with LD Score Regression
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

## Reference data - downloads

See instructions in [comorment/mixer](https://github.com/comorment/mixer/blob/main/README.md) repository for how to download reference data.
All reference data described below is based on EUR ancestry, and use ``hg19`` / ``GRCh37`` genomic build.

Here is overview of the reference files:
```
1000G_EUR_Phase3_plink/1000G.EUR.QC.[1-22].bim                      # ``--bim-file`` argument
1000G_EUR_Phase3_plink/1000G.EUR.QC.[1-22].run4.ld                  # ``--ld-file`` argument
```

Additional files for cross-trait MiXeR analysis:
```
1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep17.snps   # ``--extract`` flag
``` 

Additional files for GSA-MiXeR analysis:
```
1000G_EUR_Phase3_plink/baseline_v2.2_1000G.EUR.QC.[1-22].annot.gz   # ``--annot-file`` / ``--annot-file-test`` arguments
1000G_EUR_Phase3_plink/1000G.EUR.QC.@.[bed/bim,fam]                 # reference for MAGMA analysis, merged across chromosomes
gsa-mixer-baseline-annot_10mar2023.csv           # ``--go-file`` (baseline model)
gsa-mixer-gene-annot_10mar2023.csv               # ``--go-file`` (model)
gsa-mixer-geneset-annot_10mar2023.csv            # ``--go-file-test`` (only gene-sets)
gsa-mixer-genesetLOO-annot_10mar2023.csv         # ``--go-file-test`  (only gene-sets, with leave-one-gene-out)
gsa-mixer-hybrid-annot_10mar2023.csv             # ``--go-file-test`  (genes and gene-sets)
gsa-mixer-hybridLOO-annot_10mar2023.csv          # ``--go-file-test`  (genes and gene-sets, with leave-one-gene-out)
magma-gene-annot_10mar2023.csv                   # gsa-mixer-gene-annot_10mar2023.csv converted to MAGMA format
magma-geneset-annot_10mar2023.csv                # gsa-mixer-geneset-annot_10mar2023.csv converted to MAGMA format
```

Functional annotations are derived from [these files](https://github.com/comorment/mixer/tree/main/reference/ldsc/baselineLD_v2.2_bedfiles_only_binary) using scripts from [here](https://github.com/ofrei/eas_partitioned_ldscore). There is no need to compute LD-scores for these annotations, because MiXeR does this internally using sparse LD matrix stored in ``--ld-file`` it receives as an argument. Gene- and gene-set definitions are derived using scripts from [here](https://github.com/ofrei/genesets/blob/main/prepare_genesets_v2.ipynb).

There is an additional post-processing step needed for GSA-MiXeR analysis, producing ``.bin`` files as shown below. 
After this step loading reference is possible with ``--loadlib-file``, providing considerable speedup over passing ``--ld-file`` argument.

The following example produces ``.bin`` files for ``plsa`` analysis, yielding its own ``.bin`` file for each chromosome. The ``--savelib-file`` argument must include ``@`` symbol which will be replaced with an actual chromosome label.
In order to reduce peak memory use the command is executed through a ``for`` loop, i.e. separately for each chromosome, but it's also ok to remove ``--chr2use $chr`` option and run just a single `${MIXER_PY} plsa --savelib-file` command.
```
for chri in {1..22}; do ${MIXER_PY} plsa \
      --bim-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --use-complete-tag-indices --exclude-ranges [] --chr2use $chri \
      --savelib-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin \
      --out 1000G_EUR_Phase3_plink/1000G.EUR.QC.@; done
```

### Generate your own LD reference (cross-trait & GSA)

MiXeR reference files can be prepared from plink bfile using ``mixer.py ld`` command. It's important that the reference genotypes contain unrelated individuals only, constrained to a single population.
Note that ``@`` symbol must remain as it is in all commands, i.e. you don't need to exchange it with a specific chromosome label.

Compute LD matrices (one per chromosome), later to be used with ``--ld-file`` argument.
```
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

For full command-line reference, see ``mixer.py ld --help``.

Analyses in [GSA-MiXeR publication](https://www.nature.com/articles/s41588-024-01771-1) are based on UKB and HRC reference panen, partly shared here:
* [UKB reference](https://github.com/comorment/mixer/tree/main/reference/ukb_EUR_qc) and 
* [HRC reference](https://github.com/comorment/mixer/tree/main/reference/hrc_EUR_qc).

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
To re-generate LD matrices you will need to obtain access to individual-level data of UKB or HRC subjects, re-run QC as defined in ``about_UKB_qc.txt`` / ``about_HRC_qc.txt`` steps, and apply ``mixer.py ld``.

## FAQ

### AIC BIC interpretation

``.csv`` files generated by ``python precimed/mixer_figures.py`` commands contain AIC ([Akaike Information Criterion](https://en.wikipedia.org/wiki/Akaike_information_criterion)) and BIC ([Bayesian Information Criterion](https://en.wikipedia.org/wiki/Bayesian_information_criterion)) values. To generate AIC / BIC values you should point ``mixer_figures.py`` to json files produced by ``fit1`` or ``fit2`` steps (not those from ``test1`` or ``test2`` steps).

The idea of model selection criteia (both AIC and BIC) is to find whether the input data (in our case the GWAS summary statistics) have enough statistical power to warrant a more complex model - i.e. a model with additional free parameters that need to be optimized from the data. For example, the LDSR model has two free parameters - the slope and an intercept. The univariate MiXeR model has three parameters (``pi``, ``sig2_beta`` and ``sig2_zero``). Naturally, having an additional free parameters allows MiXeR to fit the GWAS data better compared to LDSR model, however it needs to be *substantially* better to justify an additional complexity. AIC and BIC formalize this trade-off between model's complexity and model's ability to describe input data. 

The difference between AIC and BIC is that BIC is a more conservative model selection criterion. Based on our internal use of MiXeR, it is OK discuss the resulsult if only AIC supports MiXeR model, but it's important point as a limitation that the input signal (GWAS) has borderline power to fit MiXeR model. 

For the univariate model, AIC / BIC values are described in the cross-trait MiXeR paper ([ref](https://www.nature.com/articles/s41467-019-10310-0)). A negative AIC value means that there is not enough power in the input data to justify MiXeR model as compared to LDSR model, and we do do not recommend applying MiXeR in this situation.

For the bivariate model the resulting table contains two AIC, and two BIC values, named as follows:
``best_vs_min_AIC``, ``best_vs_min_BIC`` , ``best_vs_max_AIC``, and ``best_vs_max_BIC``. They are explained below - but you may need to develop some intuition to interpret these numbers. Consider taking a look at the figure shown below, containing negative log-likelihood plots. Similar plots were presented in [ref](https://www.nature.com/articles/s41467-019-10310-0), supplementary figure 19. We use such likelihood-cost plots to visualise the performance of the ``best`` model vs ``min`` and ``max``.

First, let's interpret ``best_vs_max_AIC`` value. It uses AIC to compare two models: the ``best`` model with polygenic overlap that is shown in the venn diagram (i.e. fitted by MiXeR), versus the ``max`` model with maximum possible polygenic overlap given trait's genetic architecture (that is, in the ``max`` model the causal variants of the least  polygenic trait form a subset of the causal variants in the most polygenic trait). A positive value of ``best_vs_max_AIC`` means that ``best`` model explains the observed GWAS signal better than ``max`` model, despite its additional complexity (here additional complexity comes from the fact that MiXeR has to find an actual value of the polygenic overlap, i.e.  estimate the size of the grey area on the venn diagram).

Similarly, ``best_vs_min_AIC`` compares the ``best`` model versus model with minimal possible polygenic overlap. It is tempting to say that minimal possible polygenic overlap means the same as no shared causal variants, i.e. a case where circles on a venn diagram do not overlap. However, there is a subtle detail in our definition of the minimal possible polygenic overlap in cases where traits show some genetic correlation. Under the assumptions of cross-trait MiXeR model, the presence of genetic correlation imply non-empty polygenic overlap. In our AIC / BIC analysis, we constrain the ``best`` model to a specific value of the genetic correlation, obtained by the same procedure as in LDSR model (i.e. assuming an infinitesimal model, as described in [ref](https://www.nature.com/articles/s41467-019-10310-0). In this setting, minimal possible polygenic overlap corresponds to ``pi12 = rg * sqrt(pi1u * pi2u)``. ``best_vs_min_AIC`` is an important metric - a positive value indicates that the data shows support for existense of the polygenic overlap, beyond the minimal level need to explain observed genetic correlation between traits.

Finally, ``best_vs_max_BIC`` and ``best_vs_max_BIC`` have the same meaning as the values explained above, but using more stringent ``BIC`` criterion instead of ``AIC``.

The last panel on the following figure shows negative log-likelihood plot as a function of polygenic overlap. For the general information about maximum likelihood optimization see [here](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation). In our case we plot negative log-likelihood, that's we search for the minimum on the curve (this sometimes is refered to as a cost function, with objective to find parameters that yield the lowest cost). 

On the negative log-likelihood plot, the ``min`` model is represented by the point furthest to the left, ``max`` furthest to the right and ``best`` is the lowest point of the curve. (off note, in some other cases log-likelihood plot can be very noisy - then ``best`` is still the lowest point, but it doesn't make practical sence as it is just a very noisy estimate - and this is exactly what we want to clarify with AIC / BIC).

The minimum is obtained at ``n=1.3K`` shared causal variants - that's our ``best`` model, scores at about 25 points (lower is better). A model with least possible overlap has ``n=0`` shared causal variants - that our ``min`` model, scored at 33 points. Finally, a model with largest possible overlap has ``n=4K`` shared causal variants - that our ``max`` model, scores at `50` points. We use AIC (and BIC) to compare ``best`` model versus the other models.

![GIANT_HEIGHT_2018_UKB_vs_PGC_BIP_2016 json](https://user-images.githubusercontent.com/13171449/83339454-469edb00-a2ce-11ea-9e69-99270d94689f.png)

### Upgrade nodes from MiXeR v1.2 to v1.3

* The source code is nearly identical, but I've changed the procedure to run MiXeR which is described in this README file. Still, you need to updated the code by ``git pull``. If you updated MiXeR in early summer 2020 there will be no changes in the native C++ code (hense no need to re-compile the ``libbgmg.so`` binary), but to be safe it's best to execute ``cd <MIXER_ROOT>/src/build && cmake .. && make`` command (after ``module load`` relevant modules, as described in [this](https://github.com/precimed/mixer/blob/master/README.md#build-from-source---linux) section.)
* If you previously downloaded 1kG reference, you need to download 20 new files called ``1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.repNN.snps``. Otherwise you need to generate them with ``mixer.py snps`` command as described above.
* If you previously generated your input summary statistics constraining them to HapMap3 SNPs, you need to re-generated them without constraining to HapMap3.
* Assuming you have a SLURM script for each MiXeR analysis, you need to turn this script into a job array (``#SBATCH --array=1-20``) which executes 20 times, and use ``--extract 1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps`` flag in ``mixer.py fit1`` and ``mixer.py fit2``, also adjusting input/output file names accordingly. Example scripts are available in [scripts](https://github.com/precimed/mixer/tree/master/scripts) folder.
* To process ``.json`` files produced by MiXeR, you now first need to run a ``mixer_figures.py combine`` step as shown above. This will combine 20 ``.json`` files by averaging individual parameter estimates, and calculating standard errors.
* When you generate figures with ``mixer_figures.py one`` and ``mixer_figures.py two`` commands and use the "combined" ``.json`` files, you'll need to add ``--statistic mean std`` to your commands.
* With MiXeR v1.3 you should expect a slightly lower polygenicity estimate, mainly because of ``--maf 0.05`` constraint on the frequency of genetic variants used in the fit procedure. The rationale for ``--maf 0.05`` filter is to still have a robust selection of SNPs in the fit procedure despite not having HapMap3 filter.
* With MiXeR v1.3 you should expect a 10 fold increase in CPU resources needed per ran (20 times more runs, but each run is ~600K SNPs which is half the size of HapMap3).

### Installing Docker, Singularity, oras on your local machine

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

## Changelog

* ``v2.2.1`` - release including GSA-MiXeR functionality
* ``v1.3`` - release for Python (change HapMap3 to twenty random sets of ~600K SNPs each)
* ``v1.2`` - internal release for Python (Matlab support removed)
* ``v1.1`` - internal release for Matlab and Python
* ``v1.0`` - release for Matlab

For more detailed changelog, see [gsa-mixer/CHANGELOG.md](https://github.com/precimed/gsa-mixer/blob/main/CHANGELOG.md).
