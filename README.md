## User documentation for MiXeR analyses (univariate, bivariate and GSA-MiXeR)

> [!WARNING]
> As of Sep 11, 2025 this documentation is under active re-constructions, which will be finalized within a day or so.

This repository (https://github.com/precimed/mixer) provides user documentation for MiXeR analyses (univariate, bivariate and GSA-MiXeR).
You should use this repository for most up to date instructions on how to install and run MiXeR.
Another relevant repository is https://github.com/comorment/mixer where you will find reference data.
These two repositories contain everything you need to run MiXeR.

If you previously ran univariate or cross-trait MiXeR using ``mixer.sif`` container obtained from https://github.com/comorment/mixer project, the only step you need to do to upgrade to the latest version is to download ``gsa-mixer.sif`` container as described [here](#install), and use it instead of ``mixer.sif``. Everything else (input data, scripts performing the analysis, and output formats) remain the same, and the results should be nearly identical, except for minor bug fixes.

If you ancounter an issue, please submit a ticket: https://github.com/precimed/mixer/issues/new . In the past my response to new tickets was very bad. If you've already submitted a ticket but didn't get a response please give me a second chance to address your issue - just push a comment and tag me @ofrei , if your question is still relevant .

MiXeR source code is in https://github.com/precimed/gsa-mixer repository. This also covers univariate and cross-trait analyses, i.e. not just GSA MiXeR.
Either way, it is only relevant to developers interested to contribute pull requests to mixer code.

Additional instructions for users at the NORMENT centre are available in https://github.com/precimed/mixer_at_norment.
If the link doesn't work please reach me out by e-mail to get access.

Kind regards,
Oleksandr Frei.

## Introduction

**GSA-MiXeR** is a new technique for competitive gene-set analysis, which fits a model for gene-set heritability enrichments for complex human traits, thus allowing the quantification of partitioned heritability and fold enrichment for small gene-sets.

**Cross-trait MiXeR** is a statistical tool which quantifies polygenic overlap between complex traits irrespective of genetic correlation, using GWAS summary statistics. MiXeR results are presented as a Venn diagram of unique and shared polygenic components across traits. Cross-trait (bivariate) analysis builds on univariate analysis quantifying polygenicity of the traits.

Please cite relevant publications if you use MiXeR software in your research work.

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

To test your installation, try
```
${MIXER_PY} --help
```

The usage of ``${MIXER_PY}`` should be the same regardless of whether you use Docker or singularity version,
however for most users we recommend running through singularity container (mainly because singularity is more commonly available in HPC clusters).
If you use docker version, you may need to customize ``$DOCKER_RUN`` variable to your environment, e.g. 
you may also try replacing ``$PWD`` with ``pwd`` (same as in the above $MIXER_PY command using singularity container),
so that current working directory is mounted to the container even if you change it after defining ``MIXER_PY`` variable.

The above containers are only generated for CPUs with x86 architecture (e.g. intel or AMD CPUs), and do not support ARM architectures (for example the are not compatible with newer Macbook laptops with M1/M2/M3 chips).

## Hello-world example (cross-trait)

See [usecases/mixer_simu.md](usecases/mixer_simu.md)

## Hello-world example (GSA)

See [GSA_MIXER_README.md#getting-started-example](GSA_MIXER_README.md#getting-started-example)

## Download full reference data

See instructions in [comorment/mixer](https://github.com/comorment/mixer/blob/main/README.md) repository

## Tutorials

See [usecases/mixer_real.md](usecases/mixer_real.md) for tutorial on cross-trait MiXeR analysis
More information is in [the legacy README file](CROSS_TRAIT_README.md) describing cross-trait analysis.

See [scripts/GSA_MIXER.job](scripts/GSA_MIXER.job) for an example of script running GSA-MiXeR analysis.
More information is in [the legacy README file](GSA_MIXER_README.md) describing GSA-MiXeR analysis.

## Input data formats (cross-trait & GSA)

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

## Generate your own LD reference (cross-trait & GSA)

See [GSA_MIXER_README.md#generate-ld-reference](GSA_MIXER_README.md#generate-ld-reference). This applies to both GSA-MiXeR and cross-trait MiXeR.

## "Knowledge hub" - what you need to know about git, git-lfs, docker & signularity, oras

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

## Legacy information 

This sections (TBD) provides notes  on version changes,  notes for users upgrading from one version to another, benchmarks results comparing versions, etc;  brief history (github repos involved, reference data, approach to releasing the software versioning)

## FAQ

TBD - add domain-specific clarifications on how to interpret results, pointers to papers, etc.

