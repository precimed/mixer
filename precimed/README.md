## Contents

* [Introduction](#introduction)
* [Install MiXeR](#install-mixer)
* [Data downloads](#data-downloads)
* [Data preparation](#data-preparation)
* [Run MiXeR](#run-mixer)
* [MiXeR options](#mixer-options)
* [Visualize MiXeR results](#visualize-mixer-results)

## Introduction

This folder contains a Python port of MiXeR, wrapping the same C/C++ core as we previously used from MATLAB.
This is work in progress, but eventually it should singificantly improve user experience,
as Python allows much simpler  installation procedures,
makes it less error prone, allows to implement well-documented command-line interface (``python mixer.py --help``),
and it with visualization.

Input data for MiXeR consists of summary statistics from a GWAS, and a reference panel.
MiXeR format for summary statistics is compatible with LD Score Regression
(i.e. the ``sumstats.gz`` files), and for those users who are already familiar with ``munge_sumstats.py``
script we recommend to use LD Score Regression pipeline to prepare summary statistics.
At the same time, we encourage everyone to take a look [our own pipeline](https://github.com/precimed/python_convert/)
for processing summary statistics. For the reference panel we recommend to use 1000 Genomes Phase3 data,
pre-processed according to LD Score Regression pipeline, and available for download from LDSC website.
Further details are given in [Data downloads](#data-downloads) and [Data preparation](#data-preparation) sections.

Once you have all input data in MiXeR-compatible format you may proceed with running univariate or
and cross-trait analysis, as implemented in ``mixer.py`` command-line interface.
The results will be saved as ``.json`` files.
To visualize the results we provide a script in python, but we encourage users to write their own scripts
that process the results. Further details are given in [Run MiXeR](#run-mixer),
[MiXeR options](#mixer-options) and [Visualize MiXeR results](#visualize-mixer-results) sections.

If you encounter an issue, or have further questions, please create a
[new issue ticket](https://github.com/precimed/mixer/issues/new).

If you use MiXeR software for your research publication, please cite the following paper(s):

* for univariate analysis: D. Holland et al., Beyond SNP Heritability: Polygenicity and Discoverability Estimated for Multiple Phenotypes with a Univariate Gaussian Mixture Model, bioXriv, doi: https://doi.org/10.1101/133132
* for cross-trait analysis: O.Frei et al., Bivariate causal mixture model quantifies polygenic overlap between complex traits beyond genetic correlation, Nature Communications, 2019, https://www.nature.com/articles/s41467-019-10310-0

The MiXeR software may not be used for commercial purpose or in medical applications.
We encourage all users to familiarize themselves with US patent https://www.google.no/patents/US20150356243 "Systems and methods for identifying polymorphisms".

## Install MiXeR

### Prerequisites

* Linux environment (tested on CentOS release 6.9, Ubuntu 18.04).
* Python 3 (tested with Python 3.7), with numpy, scipy>=1.2.1, numdifftools, matplotlib_venn
* Reasonably new GCC compiler (tested with gcc/6.1.0)

If you are an experienced C++ programmer it shouldn't be difficult to compile MiXeR core in MAC or Windows.
If you've made it work, please share your changes via pull request.

### Hardware requirements

MiXeR software is very CPU and memory intensive. 
Minimal memory requirement is to have 61.5 GB of RAM available to MiXeR.
MiXeR efficiently uses multiple CPUs. We recommend to run MiXeR on a system with 16 physical cores.
When use MiXeR on a cluster, we recommend to assign the whole node to each MiXeR run.

### Install on Linux using pre-built binaries

Not available yet.

### Build from source - Linux

* Download and compile Boost libraries
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
  cmake .. -DBOOST_ROOT=$HOME/boost_1_69_0 && make bgmg -j16 
  ```

* Recommended combination of modules on different HPC systems:
  ```
  module load cmake/3.7.1 python3/3.7.0.gnu   # UiO Abel
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

* Reference panel (alternatively, download [this](https://1drv.ms/u/s!Ai1YZmdFa9ati40Inztrv_4erqcdWw?e=ixWDUe) or take it from NIRD (``/projects/NS9114K/users/oleksanf/public_file_share/mixer``). NB! Download size is ``24 GB``.
  * Run ``plink`` to calculate allele frequencies and pairwise LD r2 for each chromosome
    ```
    plink \
       --freq --threads 1 --memory 1024 \
       --bfile LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label> \
       --out LDSR/1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.<chr_label>
    plink --r2 gz --ld-window-kb 1000000 --ld-window 50000 --ld-window-r2 0.05 --threads 24 \
       --bfile LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label> \
       --out LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p05_SNPwind50k
    ```
    This step takes about 1 hour (here and below all times are measured on a server with 24 physical cores).
   * Run ``python mixer.py ld`` to convert plink output into a binary format. The following command must be run once for each chromosome. 
    Note that ``--bim`` argument is the same, e.g. ``1000G.EUR.QC.@.bim`` (with ``@``) regardless of the actual chromosome that you use in ``--plink-ld`` and ``--out``.
    ```
    python3 <MIXER_ROOT>/python/precimed/mixer.py ld \
       --lib <MIXER_ROOT>/src/build/lib/bgmglib.so
       --bim LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
       --plink-ld LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p05_SNPwind50k.ld.gz \
       --out LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p05_SNPwind50k.ld
    ```
    The conversion is single-threaded, and takes about 10 minutes for the longest chromosome.
    The output is written to ``LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p05_SNPwind50k.ld.bin`` file,
    and log details into ``LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p05_SNPwind50k.ld.log`` file:
    ```
    20181213 02:18:20.854727         Create new context (id=0)
    20181213 02:18:20.854751        >init(bim_file=LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim, frq_file=, chr_labels=, trait1_file=, trait2_file=, exclude=, extract=);
    20181213 02:18:20.857294         Construct reference from 22 files...
    20181213 02:18:22.110188         Found 141123 variants in LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.22.bim
    ...
    20181213 02:18:25.890214         Found 9997231 variants in total.
    20181213 02:18:48.881013         set_tag_indices(num_snp=9997231, num_tag=9997231);
    20181213 02:18:48.939690        >set_chrnumvec(9997231);
    20181213 02:18:48.953602        <set_chrnumvec(9997231);
    20181213 02:18:48.953646        <init(bim_file=LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim, frq_file=, chr_labels=, trait1_file=, trait2_file=, exclude=, extract=);  elapsed time 28098ms
    20181213 02:18:48.953918         Reading LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.21.p05_SNPwind50k.ld.gz...
    20181213 02:18:49.433315         Processed 100000 lines
    ...
    20181213 02:20:01.654392         Parsed 18495973 r2 values from LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.21.p05_SNPwind50k.ld.gz
    20181213 02:20:01.657704         PlinkLdFile::save_as_binary(filename=LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.21.p05_SNPwind50k.ld.bin), writing 18495973 elements...
    ```
  * Save the list of dbSNP rs# into a separate file called ``w_hm3.justrs``:
    ```
    cut -f1 w_hm3.snplist | tail -n +2 > w_hm3.justrs
    ```

## Run MiXeR

### Univariate analysis

Fit the model:
```
python3 <MIXER_ROOT>/python/precimed/mixer.py fit \
      --trait1-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --out SSGAC_EDU_2018_no23andMe_noMHC.fit \
      --extract LDSR/w_hm3.justrs --ci-alpha 0.05 \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --frq-file LDSR/1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.@.frq \
      --plink-ld-bin0 LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin \
      --lib  <MIXER_ROOT>/src/build/lib/bgmglib.so \
```

Apply the model to the entire set of SNPs, without constraining to ``LDSR/w_hm3.justrs``:
```
python3 <MIXER_ROOT>/python/precimed/mixer.py fit \
      --trait1-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --load-params-file SSGAC_EDU_2018_no23andMe_noMHC.fit.json \
      --out SSGAC_EDU_2018_no23andMe_noMHC.test \
      --fit-sequence load inflation --power-curve --qq-plots --kmax 100 \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --frq-file LDSR/1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.@.frq \
      --plink-ld-bin0 LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin \
      --lib  <MIXER_ROOT>/src/build/lib/bgmglib.so \
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
python3 <MIXER_ROOT>/python/mixer.py fit \
      --trait1-file PGC_SCZ_2014_EUR_qc_noMHC.csv.gz \
      --trait2-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --trait1-params-file PGC_SCZ_2014_EUR_qc_noMHC.fit.json \
      --trait2-params-file SSGAC_EDU_2018_no23andMe_noMHC.fit.json \
      --out PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.fit \
      --extract LDSR/w_hm3.justrs --ci-alpha 0.05 \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --frq-file LDSR/1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.@.frq \
      --plink-ld-bin0 LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin \
      --lib  <MIXER_ROOT>/src/build/lib/bgmglib.so \
```

Apply the model to the entire set of SNPs, without constraining to ``LDSR/w_hm3.justrs``:
```
python3 <MIXER_ROOT>/python/mixer.py fit \
      --trait1-file PGC_SCZ_2014_EUR_qc_noMHC.csv.gz \
      --trait2-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --load-params-file PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.fit.json \
      --out PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.test \
      --fit-sequence load inflation --power-curve --qq-plots --kmax 100 \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --frq-file LDSR/1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.@.frq \
      --plink-ld-bin0 LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin \
      --lib  <MIXER_ROOT>/src/build/lib/bgmglib.so \
```

Note that these parameters point to the results of univariate analysis for both traits, so those must be generated first.
The results will be saved ``<out_file>.json`` file.

To visualize the results:
```
python precimed/mixer_figures.py two --json <out_file>.json --out <out_file>
```

## MiXeR options

```
>python3 mixer.py fit --help

usage: mixer.py fit [-h] [--argsfile ARGSFILE] [--out OUT] [--lib LIB]
                    [--log LOG] [--bim-file BIM_FILE] [--frq-file FRQ_FILE]
                    [--plink-ld-bin PLINK_LD_BIN]
                    [--plink-ld-bin0 PLINK_LD_BIN0] [--chr2use CHR2USE]
                    [--trait1-file TRAIT1_FILE] [--trait2-file TRAIT2_FILE]
                    [--fit-sequence {load,inflation,infinitesimal,diffevo,diffevo-fast,neldermead,neldermead-fast,brute1,brute1-fast,brent1,brent1-fast} [{load,inflation,infinitesimal,diffevo,diffevo-fast,neldermead,neldermead-fast,brute1,brute1-fast,brent1,brent1-fast} ...]]
                    [--preliminary] [--extract EXTRACT] [--exclude EXCLUDE]
                    [--randprune-n RANDPRUNE_N] [--randprune-r2 RANDPRUNE_R2]
                    [--kmax KMAX] [--seed SEED] [--cache-tag-r2sum]
                    [--max-causals MAX_CAUSALS] [--r2min R2MIN]
                    [--ci-alpha CI_ALPHA] [--ci-samples CI_SAMPLES]
                    [--threads THREADS] [--tol-x TOL_X] [--tol-func TOL_FUNC]
                    [--cubature-rel-error CUBATURE_REL_ERROR]
                    [--cubature-max-evals CUBATURE_MAX_EVALS]
                    [--load-params-file LOAD_PARAMS_FILE]
                    [--trait1-params-file TRAIT1_PARAMS_FILE]
                    [--trait2-params-file TRAIT2_PARAMS_FILE]

optional arguments:
  -h, --help            show this help message and exit
  --argsfile ARGSFILE   file with additional command-line arguments
  --out OUT             prefix for the output files
  --lib LIB             path to libbgmg.so plugin
  --log LOG             file to output log, defaults to <out>.log
  --bim-file BIM_FILE   Plink bim file. Defines the reference set of SNPs used
                        for the analysis. Marker names must not have
                        duplicated entries. May contain simbol '@', which will
                        be replaced with the actual chromosome label.
  --frq-file FRQ_FILE   Plink frq file (alleles frequencies). May contain
                        simbol '@', similarly to --bim-file argument.
  --plink-ld-bin PLINK_LD_BIN
                        File with linkage disequilibrium information,
                        converted from plink format as described in the
                        README.md file. May contain simbol '@', similarly to
                        --bim-file argument.
  --plink-ld-bin0 PLINK_LD_BIN0
                        File with linkage disequilibrium information in an old
                        format (deprecated)
  --chr2use CHR2USE     Chromosome ids to use (e.g. 1,2,3 or 1-4,12,16-20).
                        Chromosome must be labeled by integer, i.e. X and Y
                        are not acceptable.
  --trait1-file TRAIT1_FILE
                        GWAS summary statistics for the first trait.
  --trait2-file TRAIT2_FILE
                        GWAS summary statistics for the first trait.
                        Specifying this argument triggers cross-trait
                        analysis.
  --fit-sequence {load,inflation,infinitesimal,diffevo,diffevo-fast,neldermead,neldermead-fast,brute1,brute1-fast,brent1,brent1-fast} [{load,inflation,infinitesimal,diffevo,diffevo-fast,neldermead,neldermead-fast,brute1,brute1-fast,brent1,brent1-fast} ...]
                        Specify fit sequence: 'load' reads previosly fitted
                        parameters from a file (--load-params-file); 'diffevo'
                        performs a single iteration of differential evolution,
                        which is the right way to setup an initial
                        approximation; 'neldermead' applies Nelder-Mead
                        downhill simplex search; 'brute1' applies only to
                        bivariate optimization; it performs brute-force for
                        one-dimentional optimization, searching optimal pi12
                        value constrained on genetic correlation (rg) and
                        intercept (rho0); 'brent1' is similar to brute1, but
                        it uses brent method (inverse parabolic
                        interpolation); 'inflation' fits sig2zero (univariate)
                        and rho_zero (bivariate), using fast cost function;
                        this is quite special optimization step, typically
                        useful for adjusting inflation parameters to another
                        reference; 'infinitesimal' fits a model with pi1=1
                        (univariate) or pi12=1 (bivariate) constrains, using
                        fast cost function; this is quite special optimization
                        step, typically used internally for AIC/BIC
                        computation; Note that bivariate fit is always
                        constrained on univariate parameters, except for
                        'inflation' fit which adjust rho_zero and sig2_zero.
                        The '...-fast' optimizations use fast cost function.
                        Note that univariate optimization uses 'convolve' cost
                        calculator, bivariate optimization uses 'sampling'
                        cost calculator. Typical univariate sequence:
                        'diffevo-fast neldermead'Typical bivariate sequence:
                        'diffevo neldermead brute1 brent1'
  --preliminary         perform an additional run using fast model with
                        'diffevo-fast nelderead-fast' to generate preliminary
                        data. After preliminary run fit sequence is applied
                        from scratch using full model.
  --extract EXTRACT     File with variants to include in the fit procedure
  --exclude EXCLUDE     File with variants to exclude from the fit procedure
  --randprune-n RANDPRUNE_N
                        Number of random pruning iterations
  --randprune-r2 RANDPRUNE_R2
                        Threshold for random pruning
  --kmax KMAX           Number of sampling iterations
  --seed SEED           Random seed
  --cache-tag-r2sum     enable tag-r2sum caching
  --max-causals MAX_CAUSALS
                        upper limit for the total number of causal variants in
                        the reference; a number between 0 and 1 represents a
                        fraction of the total number SNPs in the reference
  --r2min R2MIN         r2 values below this threshold will contribute via
                        infinitesimal model
  --ci-alpha CI_ALPHA   significance level for the confidence interval
                        estimation
  --ci-samples CI_SAMPLES
                        number of samples in uncertainty estimation
  --threads THREADS     specify how many threads to use (concurrency). None
                        will default to the total number of CPU cores.
  --tol-x TOL_X         tolerance for the stop criteria in fminsearch
                        optimization.
  --tol-func TOL_FUNC   tolerance for the stop criteria in fminsearch
                        optimization.
  --cubature-rel-error CUBATURE_REL_ERROR
                        relative error for cubature stop criteria (applies to
                        'convolve' cost calculator).
  --cubature-max-evals CUBATURE_MAX_EVALS
                        max evaluations for cubature stop criteria (applies to
                        'convolve' cost calculator). Bivariate cubature
                        require in the order of 10^4 evaluations and thus is
                        much slower than sampling, therefore it is not exposed
                        via mixer.py command-line interface.
  --load-params-file LOAD_PARAMS_FILE
                        initial params for the optimization.
  --trait1-params-file TRAIT1_PARAMS_FILE
                        univariate params for the first trait (for the cross-
                        trait analysis only).
  --trait2-params-file TRAIT2_PARAMS_FILE
                        univariate params for the second trait (for the cross-
                        trait analysis only).
```

### Memory usage
 
TBD. MiXeR is still using a lot of memory, but we are working on making it better.

## Visualize MiXeR results

TBD.

### MiXeR results format

TBD.

All MiXeR results are stored in a single ``.json`` file.
