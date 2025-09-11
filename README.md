## User documentation for MiXeR analyses (univariate, bivariate and GSA-MiXeR)

> [!WARNING]
> As of Sep 11, 2025 this documentation is under active re-constructions, which will be finalized within a day or so.

This repository (https://github.com/precimed/mixer) provides user documentation for MiXeR analyses (univariate, bivariate and GSA-MiXeR).
You should use this repository for most up to date instructions on how to install and run MiXeR analyses.
Another relevant repository is https://github.com/comorment/mixer where you will find reference data.
These two repositories contain everything you need to run univariate, bivariate or GSA-MiXeR analyses.

If you ancounter an issue, please submit a ticket: https://github.com/precimed/mixer/issues/new . In the past my response to new tickets was very bad. If you've already submitted a ticket but didn't get a response please give me a second chance to address your issue - just push a comment and tag me @ofrei , if your question is still relevant .

Source code for MiXeR tools is in https://github.com/precimed/gsa-mixer repository.
It is only relevant to developers interested to contribute pull requests to mixer code.

Additional instructions for users at the NORMENT centre are available in https://github.com/precimed/mixer_at_norment
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

## Install (cross-trait & GSA)

See [GSA_MIXER_README.md#install-gsa-mixer](GSA_MIXER_README.md#install-gsa-mixer)

## Hello-world example (cross-trait)

See [usecases/mixer_simu.md](usecases/mixer_simu.md)

## Hello-world example (GSA)

See [GSA_MIXER_README.md#getting-started-example](GSA_MIXER_README.md#getting-started-example)

## Download full reference data

See instructions in [comorment/mixer](https://github.com/comorment/mixer/blob/main/README.md) repository

## Tutorials

See [usecases/mixer_real.md](usecases/mixer_real.md) for tutorial on cross-trait MiXeR analysis
More information is in [the legacy README file](CROSS_TRAIT_README.md) describing cross-trait analysis.

See [v](scripts/GSA_MIXER.job) for an example of script running GSA-MiXeR analysis.
More information is in [the legacy README file](GSA_MIXER_README.md) describing GSA-MiXeR analysis.

## Input data formats (cross-trait & GSA)

See [GSA_MIXER_README.md#input-data-formats](GSA_MIXER_README.md#input-data-formats). This applies to both GSA-MiXeR and cross-trait MiXeR.

## Generate your own LD reference (cross-trait & GSA)

See [GSA_MIXER_README.md#generate-ld-reference](GSA_MIXER_README.md#generate-ld-reference). This applies to both GSA-MiXeR and cross-trait MiXeR.

## "Knowledge hub" - what you need to know about git, git-lfs, docker & signularity, oras

TBD

## Legacy information 

This sections (TBD) provides notes  on version changes,  notes for users upgrading from one version to another, benchmarks results comparing versions, etc;  brief history (github repos involved, reference data, approach to releasing the software versioning)

## FAQ

TBD - add domain-specific clarifications on how to interpret results, pointers to papers, etc.

