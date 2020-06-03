## Contents

* [Introduction](#introduction)
* [Install MiXeR](#install-mixer)
* [Data downloads](#data-downloads)
* [Data preparation](#data-preparation)
* [Run MiXeR](#run-mixer)
* [MiXeR options](#mixer-options)
* [Visualize MiXeR results](#visualize-mixer-results)
* [AIC and BIC interpretation](#aic-bic-interpretation)

## Introduction

This folder contains a Python port of MiXeR, wrapping the same C/C++ core as we previously used from MATLAB.
This is work in progress, but eventually it should singificantly improve user experience,
as Python allows much simpler  installation procedures,
makes it less error prone, allows to implement well-documented command-line interface (``python mixer.py --help``),
and provide visualization.

Input data for MiXeR consists of summary statistics from a GWAS, and a reference panel.
MiXeR format for summary statistics is compatible with LD Score Regression
(i.e. the ``sumstats.gz`` files), and for those users who are already familiar with ``munge_sumstats.py``
script we recommend to use LD Score Regression pipeline to prepare summary statistics.
At the same time, we encourage everyone to take a look [our own pipeline](https://github.com/precimed/python_convert/)
for processing summary statistics. For the reference panel we recommend to use 1000 Genomes Phase3 data,
pre-processed according to LD Score Regression pipeline, and available for download from LDSC website.
Further details are given in [Data downloads](#data-downloads) and [Data preparation](#data-preparation) sections.

Once you have all input data in MiXeR-compatible format you may proceed with running univariate (``fit1``, ``test1``)
and cross-trait (``fit2``, ``test2``) analyses, as implemented in ``mixer.py`` command-line interface.
The results will be saved as ``.json`` files.
To visualize the results we provide a script in python, but we encourage users to write their own scripts
that understand the structure of ``.json`` files, process the results.
Further details are given in [Run MiXeR](#run-mixer),
[MiXeR options](#mixer-options) and [Visualize MiXeR results](#visualize-mixer-results) sections.

If you encounter an issue, or have further questions, please create a
[new issue ticket](https://github.com/precimed/mixer/issues/new).

If you use MiXeR software for your research publication, please cite the following paper(s):

* for univariate analysis: D. Holland et al., Beyond SNP Heritability: Polygenicity and Discoverability Estimated for Multiple Phenotypes with a Univariate Gaussian Mixture Model, bioXriv, doi: https://doi.org/10.1101/133132
* for cross-trait analysis: O.Frei et al., Bivariate causal mixture model quantifies polygenic overlap between complex traits beyond genetic correlation, Nature Communications, 2019, https://www.nature.com/articles/s41467-019-10310-0

The MiXeR software may not be used for commercial purpose or in medical applications.
We encourage all users to familiarize themselves with US patent https://www.google.no/patents/US20150356243 "Systems and methods for identifying polymorphisms".

MiXeR versions:

* ``v1.0`` - public release for Matlab
* ``v1.1`` - internal release for Matlab and Python
* ``v1.2`` - internal release for Python (Matlab support removed)

## Install MiXeR

### Prerequisites

* Linux environment (tested on CentOS release 6.9, Ubuntu 18.04).
* Python 3 (tested with Python 3.7), with numpy, scipy>=1.2.1, numdifftools, matplotlib_venn
* Reasonably new GCC compiler (tested with gcc/6.1.0) or an Intel C++ compiler (tested with intel-2018b)
* Boost libraries (https://www.boost.org/, tested with 1.68.0)

If you are an experienced C++ programmer it shouldn't be difficult to compile MiXeR core in MAC or Windows.
If you've made it work, please share your changes via pull request.

### Hardware requirements

MiXeR software is very CPU intensive. 
Minimal memory requirement is to have 32 GB of RAM available to MiXeR.
MiXeR efficiently uses multiple CPUs.
We recommend to run MiXeR on a system with at least 16 physical cores.

### Install on Linux using pre-built binaries

Not available yet.

### Build from source - Linux

The exact steps depend  on your build environment. 
* If you work in HPC environment with modules system, you can load some existing combination of modules that include Boost libraries:
  ```
  module load Boost/1.68.0-intel-2018b-Python-3.6.6 Python/3.6.6-intel-2018b CMake/3.12.1  # SAGA (intel)
  module load Boost/1.71.0-GCC-8.3.0 Python/3.7.4-GCCcore-8.3.0 CMake/3.12.1               # SAGA (gcc)  
  ```
* Alternatively, you may download and compile Boost libraries yourself:
  ```
  cd ~ && wget https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz 
  tar -xzvf boost_1_69_0.tar.gz && cd boost_1_69_0
  ./bootstrap.sh --with-libraries=program_options,filesystem,system,date_time
  ./b2 --clean && ./b2 --j12 -a
  ```
* Clone and compile MiXeR repository
  ```
  cd ~ && git clone --recurse-submodules -j8 https://github.com/precimed/mixer.git
  mkdir mixer/src/build && cd mixer/src/build
  cmake .. && make bgmg -j16                                   # if you use GCC compiler
  CC=icc CXX=icpc cmake .. && make bgmg -j16                   # if you use Intel compiler
  cmake .. -DBOOST_ROOT=$HOME/boost_1_69_0 && make bgmg -j16   # if you use locally compiled boost
  
  ```

## Data downloads


* Summary statistics, for example
  * Schizophrenia GWAS by Ripke at al. (2014) (search for *Download 49 EUR samples* [here](https://www.med.unc.edu/pgc/results-and-downloads))
  * Educational Attainment GWAS by Lee et al. (2018) (search for *GWAS_EA_excl23andMe.txt* [here](https://www.thessgac.org/data))

* Download reference data from [this URL](https://data.broadinstitute.org/alkesgroup/LDSCORE/)
  ```
  wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz
  wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
  tar -xzvf 1000G_Phase3_plinkfiles.tgz
  bzip2 -d w_hm3.snplist.bz2
  ```
  
## Data preparation

* Summary statistics (NIRD: ``/projects/NS9114K/MMIL/SUMSTAT/TMP/nomhc/``)
  * MiXeR recognizes summary statistics in LDSC format as described [here](https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format). In brief, each trait must be represented as a single table containing columns SNP, N, Z, A1, A2. Thus, it is possible to use ``munge_sumstats.py`` script as described [here](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability#step-1-download-the-data). This might be convenient for users who are already familiar with LDSR functionality.
  * However, we recommed to use our own scripts to pre-process summary statistcs (clone from [here](https://github.com/precimed/python_convert)):
    ```
    python sumstats.py csv --sumstats daner_PGC_SCZ49.sh2_mds10_1000G-frq_2.gz --out PGC_SCZ_2014_EUR.csv --force --auto --head 5 --ncase-val 33640 --ncontrol-val 43456 
    python sumstats.py zscore --sumstats PGC_SCZ_2014_EUR.csv | \
    python sumstats.py qc --exclude-ranges 6:26000000-34000000 --max-or 1e37 | \
    python sumstats.py neff --drop --factor 4 --out PGC_SCZ_2014_EUR_qc_noMHC.csv --force 
    gzip PGC_SCZ_2014_EUR_qc_noMHC.csv

    python sumstats.py csv --sumstats GWAS_EA_excl23andMe.txt.gz --out SSGAC_EDU_2018_no23andMe.csv --force --auto --head 5 --n-val 766345
    python sumstats.py zscore --sumstats SSGAC_EDU_2018_no23andMe.csv | \
    python sumstats.py qc --exclude-ranges 6:26000000-34000000 --out SSGAC_EDU_2018_no23andMe_noMHC.csv --force
    gzip SSGAC_EDU_2018_no23andMe_noMHC.csv  
    ```
  * We note that for case/control ``munge_sumstats.py`` generate sample size as a sum ``n = ncase + ncontrol``. We recommend to use ``neff = 4 / (1/ncase + 1/ncontrol)`` to account for imbalanced classes. Additionaly, we recommend to keep summary statistics for the entire set of SNPs available in GWAS, without filtering by HapMap3 SNPs). HapMap3 constraint can be applied later during fit procedure.

* Reference panel (alternatively, download [this](https://1drv.ms/u/s!Ai1YZmdFa9ati40Inztrv_4erqcdWw?e=ixWDUe) or take it from NIRD (``/projects/NS9114K/MMIL/SUMSTAT/LDSR/1000G_EUR_Phase3_plink``). NB! Download size is ``24 GB``.
   * Run ``python mixer.py ld`` to calculate linkage disequilibrium information in a genotype reference panel. The following command must be run once for each chromosome. 
    ```
    python3 <MIXER_ROOT>/precimed/mixer.py ld \
       --lib <MIXER_ROOT>/src/build/lib/libbgmg.so \
       --bfile LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label> \
       --out LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.run4.ld \
       --r2min 0.05 --ldscore-r2min 0.05 --ld-window-kb 30000
    ```
    The output is written to ``LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.run4.ld`` file,
    and log details into ``LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.run4.ld.log`` file.
    The resulting files contain information about alleling LD r2 correlations,
    LD scores, and allele frequencies of the variants in the reference panel passed as ``--bfile`` argument.
    The files DO NOT contain any individual-level information.
    When you store the resulting ``.ld`` file, it is important to keep it along side with corresponding ``.bim`` file,
    as information about marker name (SNP rs#), chromosome,  position, and alleles (A1/A2) is NOT encoded in ``.ld`` file.
    
  * Save the list of dbSNP rs# into a separate file called ``w_hm3.justrs``:
    ```
    cut -f1 w_hm3.snplist | tail -n +2 > w_hm3.justrs
    ```

## Run MiXeR

### Univariate analysis

Fit the model:
```
python3 <MIXER_ROOT>/precimed/mixer.py fit1 \
      --trait1-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --out SSGAC_EDU_2018_no23andMe_noMHC.fit \
      --extract LDSR/w_hm3.justrs \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  <MIXER_ROOT>/src/build/lib/libbgmg.so \
```

Apply the model to the entire set of SNPs, without constraining to ``LDSR/w_hm3.justrs``:
```
python3 <MIXER_ROOT>/precimed/mixer.py test1 \
      --trait1-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --load-params-file SSGAC_EDU_2018_no23andMe_noMHC.fit.json \
      --out SSGAC_EDU_2018_no23andMe_noMHC.test \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  <MIXER_ROOT>/src/build/lib/libbgmg.so \
```

The results will be saved ``<out_file>.json`` file.
Repeat the above analysis for the second trait (``PGC_SCZ_2014_EUR_qc_noMHC.csv.gz``).

To visualize the results:
```
python precimed/mixer_figures.py one --json <out_file>.json --out <out_file>
```


### Bivariate (cross-trait) analysis

Fit the model:
```
python3 <MIXER_ROOT>/python/mixer.py fit2 \
      --trait1-file PGC_SCZ_2014_EUR_qc_noMHC.csv.gz \
      --trait2-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --trait1-params-file PGC_SCZ_2014_EUR_qc_noMHC.fit.json \
      --trait2-params-file SSGAC_EDU_2018_no23andMe_noMHC.fit.json \
      --out PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.fit \
      --extract LDSR/w_hm3.justrs \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  <MIXER_ROOT>/src/build/lib/libbgmg.so \
```

Apply the model to the entire set of SNPs, without constraining to ``LDSR/w_hm3.justrs``:
```
python3 <MIXER_ROOT>/python/mixer.py test2 \
      --trait1-file PGC_SCZ_2014_EUR_qc_noMHC.csv.gz \
      --trait2-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --load-params-file PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.fit.json \
      --out PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.test \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  <MIXER_ROOT>/src/build/lib/libbgmg.so \
```

Note that these parameters point to the results of univariate analysis for both traits, so those must be generated first.
The results will be saved ``<out_file>.json`` file.

To visualize the results:
```
python precimed/mixer_figures.py two --json <out_file>.json --out <out_file>
```

## MiXeR options

Run ``--help`` commands to list available options and their description.
```
python3 mixer.py ld --help
python3 mixer.py fit1 --help
python3 mixer.py test1 --help
python3 mixer.py fit2 --help
python3 mixer.py test2 --help
python3 mixer.py perf --help
```

## Visualize MiXeR results

The resulting ``.json`` files can be converted to figures and ``.csv`` tables via the following commands (``one`` for univariate, ``two`` for bivariate; each of these commands accept ``.json`` files from ``fit`` and ``test`` steps).

```
python precimed/mixer_figures.py one --json <out_file>.json --out <out_file>
python precimed/mixer_figures.py two --json <out_file>.json --out <out_file>
python precimed/mixer_figures.py two --json-fit <out_file_fit>.json --json-test <out_file_test>.json --out <out_file>
```

For the ``two`` command, instead of ``--json``, it is possible to specify ``--json-fit`` and ``--json-test`` separately.
This allows to combine negative log-likelihood plot (available in fit2 only) and QQ plots (available in test2 only).
Note that all ``--json`` accept wildcards (``*``) or a list of multiple files. This allows to generate ``.csv`` tables
containing results from multiple MiXeR runs. 

### MiXeR results format

MiXeR produces the following results, as described in the original publication.

* Univariate ``.csv`` file, including model parameters and AIC/BIC values,
* Univariate QQ plots
* Univariate MAF- and LD-stratified QQ plots
* Univariate Power curves
* Bivariate ``.csv`` file, including model parameters, AIC/BIC values and Dice coefficient
* Bivariate stratified QQ plots (cross-trait enrichment)
* Bivariate density plots
* Bivariate negative log-likelihood function

TBD - provide more details.

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


