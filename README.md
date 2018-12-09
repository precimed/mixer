## Contents

* [Introduction](#introduction)
* [Install MiXeR](#install-mixer)
* [Data downloads](#data-downloads)
* [Data preparation](#data-preparation)
* [Run MiXeR](#run-mixer)
* [Visualize MiXeR results](#visualize-mixer-results)

## Introduction

Mixer code is generally implemented in Matlab, but some routines were coded in native C or C++ language to give better performance. Therefore, to run MiXeR one needs to either compile C/C++ code, or install per-built binaries which MiXeR depends on. Further details are available in [Install MiXeR](#install-mixer) section.

Input data for MiXeR consists of summary statistics from a GWAS, and a reference panel. MiXeR format for summary statistics is compatible with LD Score Regression (i.e. the ``sumstats.gz`` files), and for those users who are already familiar with ``munge_sumstats.py`` script we recommend to use LD Score Regression pipeline to prepare summary statistics. At the same time, we encourage everyone to take a look [our own pipeline](https://github.com/precimed/python_convert/) for processing summary statistics. For the reference panel we recommend to use 1000 Genomes Phase3 data, pre-processed according to LD Score Regression pipeline, and available for download from LDSC website. Further details are given in [Data downloads](#data-downloads) and [Data preparation](#data-preparation) sections.

Once you have all input data in MiXeR-compatible format you may proceed with running univariate analysis ([UGMG_cpp_run_simple.m](UGMG_cpp_run_simple.m) script) and cross-trait analysis ([BGMG_cpp_run_simple.m](BGMG_cpp_run_simple.m) script). The results will be saved as ``.json`` files. To visualize the results we provide a script in python, but we encourage users to write their own scripts that process the results. Further details are given in [Run MiXeR](#run-mixer) and [Visualize MiXeR results](#visualize-mixer-results) sections.

If you encounter an issue, or have further questions, please create a [new issue ticket](https://github.com/precimed/mixer/issues/new).

If you use MiXeR software for your research publication, please cite the following paper(s):

* for univariate analysis: D. Holland et al., Beyond SNP Heritability: Polygenicity and Discoverability Estimated for Multiple Phenotypes with a Univariate Gaussian Mixture Model, bioXriv, doi: https://doi.org/10.1101/133132
* for cross-trait analysis: O.Frei et al., Bivariate causal mixture model quantifies polygenic overlap between complex traits beyond genetic correlation, bioXriv, doi: https://doi.org/10.1101/240275 

The MiXeR software may not be used for commercial purpose or in medical applications.
We encourage all users to familiarize themselves with US patent https://www.google.no/patents/US20150356243 "Systems and methods for identifying polymorphisms".

## Install MiXeR

### Prerequisites

* MiXeR was tested on Linux (CentOS release 6.9) and Windows 10 (build 10.0.17134) operating systems
* Matlab (Release R2017a). Other versions of the Matlab may work as well, but as of today they are not guarantied to be compatible with pre-built MiXeR binaries (C/C++ code).
* Python 2.7 (for LD score regression)
* Python 3.5 (for MiXeR results visualization)
* [plink 1.9](https://www.cog-genomics.org/plink2) for LD structure estimation
* ``munge_sumstats.py`` script from [LDSC](https://github.com/bulik/ldsc) to process summary statistics

### Hardware requirements

MiXeR software is very CPU and memory intensive. 
Minimal memory requirement is to have 61.5 GB of RAM available to MiXeR.
MiXeR efficiently uses multiple CPUs. We recommend to run MiXeR on a system with 16 physical cores.
When use MiXeR on a cluster, we recommend to assign the whole node to each MiXeR run.

### Install on Linux using pre-built binaries

* Download "Linux_x64.tar.gz" file from the latest MiXeR release (https://github.com/precimed/mixer/releases)
* Extract "Linux_x64.tar.gz" to a new folder. Below we refer to this folder as ``MIXER_ROOT``.
* Test that MiXeR executable runs smoothly:
  * Start new command line
  * Change active folder to ``MIXER_ROOT``
  * Run ``bin/bgmg-cli`` command, and validate that it produces output:
    ```
    >bin/bgmg-cli --help
    BGMG v0.9.0 - Univariate and Bivariate causal mixture models for GWAS:
      -h [ --help ]         produce this help message
      ...
    ```
* Test that MiXeR C++ plugin is loaded correctly
  * open a new instance matlab
  * change active folder to ``MIXER_ROOT``
  * execute ``test_mixer_plugin`` command, and check that it did not crash
  * check that ``test_mixer_plugin`` have created a new file named ``test_mixer_plugin.bgmglib.log`` containing the following lines:
    ```
    20181208 20:16:56.281717	============= new session =============
    20181208 20:16:56.281717	=mixer plugin setup success
    ```
  * If ``test_mixer_plugin`` crashes your matlab instance, try to run the following command before starting MATLAB:
    ```
    export "LD_LIBRARY_PATH=$MIXER_ROOT/lib:$LD_LIBRARY_PATH"
    ```
	
### Install on Windows using pre-built binaries

* Download "Windows_x64.7z" file from the latest MiXeR release (https://github.com/precimed/mixer/releases)
* Extract "Windows_x64.7z" to a new folder. Below we refer to this folder as ``MIXER_ROOT``.
* Test that MiXeR executable runs smoothly, as described in [Install on Linux using pre-built binaries](#install-on-linux-using-pre-built-binaries) section
* Test that MiXeR C++ plugin is loaded correctly, as described in [Install on Linux using pre-built binaries](#install-on-linux-using-pre-built-binaries) section

### Build from source - Linux

Preliminary notes are available in [src/README.md](src/README.md).
		
### Build from source - Windows

Preliminary notes are available in [src/README.md](src/README.md).

## Data downloads

* Summary statistics, for example
  * Schizophrenia GWAS by Ripke at al. (2014) (search for *Download 49 EUR samples* [here](https://www.med.unc.edu/pgc/results-and-downloads))
  * Educational Attainment GWAS by Lee et al. (2018) (search for *GWAS_EA_excl23andMe.txt* [here](https://www.thessgac.org/data))

* Download reference data from [this URL](https://data.broadinstitute.org/alkesgroup/LDSCORE/)
  ```
  wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz
  wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
  ```

## Data preparation

* Summary statistics
  * MiXeR recognizes summary statistics in LDSC format as described [here](https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format). In brief, each trait must be represented as a single table containing columns SNP, N, Z, A1, A2. 
  * We recommend to use ``munge_sumstats.py`` script as described [here](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability#step-1-download-the-data).
  * We note that for case/control ``munge_sumstats.py`` generate sample size as a sum ``n = ncase + ncontrol``. We recommend to use ``neff = 4 / (1/ncase + 1/ncontrol)`` to account for imbalanced classes.

* Reference panel
  * Run ``plink`` to calculate allele frequencies and pairwise LD r2 for each chromosome
    ```
    plink \
       --freq --threads 1 --memory 1024 \
       --bfile LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label> \
       --out LDSR/1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.<chr_label>
    plink --r2 gz --ld-window-kb 1000000 --ld-window 50000 --ld-window-r2 0.05 --threads 24 \
       --bfile LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label> \
       --out LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p01_SNPwind50k
    ```
  * Run ``bgmg-cli`` to convert plink output into a binary format. The following command must be run once for each chromosome. 
    Note that ``--bim`` argument must be the same, i.e. point to ``1000G.EUR.QC.@.bim`` regardless of the actual chromosome that you use in ``--plink-ld`` and ``--out``.
    ```
    bin/bgmg-cli \
       --bim LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
       --plink-ld LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p01_SNPwind50k.ld.gz \
       --out LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p01_SNPwind50k.ld.bin
    ```
  * Save the list of dbSNP rs# into a separate file called ``w_hm3.justrs``:
    ```
    cut -f1 w_hm3.snplist | tail -n +2 > w_hm3.justrs
    ```

## Run MiXeR

### Univariate analysis

To start univariate analysis one should start matlab, configure parameters as listed below, and execute [UGMG_cpp_run_simple](UGMG_cpp_run_simple.m) script.
The recommended way is to call ``cd $MIXER_ROOT && matlab -nodisplay -nosplash -nodesktop -r "UGMG_params_fit; UGMG_cpp_run_simple; exit;``, where ``UGMG_params_fit.m`` is a file that looks as follows:
```
trait1_file='SSGAC_EDU_2018_no23andMe_noMHC.sumstats.gz';
out_file='SSGAC_EDU_2018_no23andMe_noMHC.fit';
bim_file='LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim';
frq_file='1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.@.frq';
plink_ld_bin='LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin'; chr_labels = 1:22;
init_from_params_file=''; extract='w_hm3.justrs';
bgmg_shared_library='lib/libbgmg.so';
bgmg_shared_library_header='bgmg_matlab.h';
kmax=20000; max_causal_fraction=0.03; z1max=nan; z2max=nan;
cache_tag_r2sum=0; r2min=0.0;
randprune_r2=0.1; randprune_n=64;
DO_FIT_UGMG=1; FIT_FULL_MODEL=1; CI_ALPHA=0.05; SEED=123;
POWER_PLOT=1; POWER_PLOT_DOWNSCALE=100;
QQ_PLOT=1; QQ_PLOT_DOWNSCALE=100;
QQ_PLOT_BINS=1; QQ_PLOT_BINS_DOWNSCALE=50;
```
The results will be saved ``<out_file>.json`` file.

The above parameters imply that the model is fitted on HapMap3 SNPs, as specified by ``extract='w_hm3.justrs'``.
To  test fitted  parameters on the entire set of SNPs one may change parameters as follows:
```
DO_FIT_UGMG=0; kmax=100; extract=''; 
init_from_params_file='SSGAC_EDU_2018_no23andMe_noMHC.fit.params.mat'
out_file='SSGAC_EDU_2018_no23andMe_noMHC.test';
```

### Bivariate (cross-trait) analysis

To start cross-trait analysis one should start matlab, configure parameters as listed below, and execute [BGMG_cpp_run_simple](BGMG_cpp_run_simple.m) script. The recommended way is to call ``cd $MIXER_ROOT && matlab -nodisplay -nosplash -nodesktop -r "BGMG_params_fit; BGMG_cpp_run_simple; exit;``, where ``BGMG_params_fit.m`` is a file that looks as follows:
```
trait1_file='PGC_SCZ_2014_EUR_qc_noMHC.sumstats.gz';
trait2_file='SSGAC_EDU_2018_no23andMe_noMHC.sumstats.gz';
out_file='PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.fit';
trait1_params_file='PGC_SCZ_2014_EUR_qc_noMHC.fit.params.mat';
trait2_params_file='SSGAC_EDU_2018_no23andMe_noMHC.fit.params.mat';
bim_file='LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim';
frq_file='LDSR/1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.@.frq';
plink_ld_bin='LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin'; chr_labels = 1:22;
init_from_params_file=''; extract='w_hm3.justrs';
bgmg_shared_library='lib/libbgmg.so';
bgmg_shared_library_header='bgmg_matlab.h';
kmax=20000; max_causal_fraction=0.02; z1max=nan; z2max=nan;
cache_tag_r2sum=0; r2min=0.0;
randprune_r2=0.1; randprune_n=64;
DO_FIT_BGMG=1; FIT_FULL_MODEL=1; CI_ALPHA=0.05; SEED=123;
STRATIFIED_QQ_PLOT=0;STRATIFIED_QQ_PLOT_DOWNSCALE=100;qq_zgrid_lim=25;qq_zgrid_step=0.05;
```
Note that these parameters point to the results of univariate analysis for both traits, so those must be generated first.
The results will be saved ``<out_file>.json`` file.

To test fitted  parameters on the entire set of SNPs one may change parameters as follows:
```
DO_FIT_BGMG=0; kmax=100; extract='';
init_from_params_file='PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.fit.params.mat';
out_file='PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.test';
```

### Sequence of operations

MiXeR performs the following sequence of operations. Some stems might be skipped, depending on user-provided options.
Operations specific to univariate analysis are marked with [UGMG] tag, 
those specific to cross-trait analysis - with [BGMG] tag, and common operations - with [BOTH] tag.

* [BOTH] Load and initialize ``bgmglib`` (native c/c++ plugin), 3rd party components (``DERIVESTsuite``), s
* [BOTH] Load reference set of SNPs (SNP/CHR/BP/A1/A2) from .bim file(s), plink format
* [BOTH] Load allele frequencies for the reference panel from .frq file(s), plink format
* [BOTH] Load trait(s) data (z-scores, and per-snp sample size) from .sumstats.gz LDSR-formatted file
* [BOTH] Load LD structure from mixer-formatted ``.ld.bin`` files (derived from plink ``.ld.gz`` using ``bgmg-cli``)
* [BOTH] Generate weights for all SNPs defined across trait(s) based on random pruning
* [BGMG] load univariate parameters for each trait, to constrain bivariate optimization
* [BOTH] Setup *initial approximation* of model parameters. Two options are available: load params from a previous run, or find initial approximation using *fast model* (gaussian approximation of the log-likelihood cost function)
* [UGMG] fit univariate parameters using full model for GWAS z scores, estimate standard errors
* [UGMG] produce QQ plots
* [UGMG] produce partitioned QQ plots (bins of MAF and LD score)
* [UGMG] produce power curves (proportion of heritability explained as a function of GWAS sample size)
* [BGMG] fit bivariate parameters using full model for GWAS z scores, estimate standard errors
* [BGMG] produce stratified QQ plots
* [BOTH] Save results to <out_file>.[json, mat, pdf, log]

### List of all MiXeR options

Options specific to univariate analysis are marked with [UGMG] tag, 
those specific to cross-trait analysis - with [BGMG] tag, 
common operations - with [BOTH] tag.

* input data files (all paths could be absolute or relative)
  * [BOTH] ``trait1_file``, required -- path to ``.sumstats.gz`` file with summary statistics in LDSR format
  * [BGMG] ``trait2_file``, required -- second trait for bivariate analysis
  * [BGMG] ``trait1_params_file``, required -- path to ``<out_tag>.params.mat`` file from an existing univariate analysis of the first trait
  * [BGMG] ``trait2_params_file``, required -- same as ``trait1_params_file`` but for the second trait
  * [BOTH] ``init_from_params_file``, optional -- path to ``<out_tag>.params.mat`` file from an existing analysis (either univariate or bivariate).
    In cross-trait analysis, when ``init_from_params_file`` is specified, the ``trait1_params_file`` and ``trait2_params_file`` params are not in use (and thus are optional).
  * [BOTH] ``bim_file``, required -- path to plink ``.bim`` files that define the reference set of SNPs. The files may be split per chromosome, 
    in which case ``bim_file`` argument must have ``@`` sign indicating the location of chromosome label.
	The list of chromosome labels must be provided via ``chr_labels`` parameter.
  * [BOTH] ``frq_file``, required -- path to plink ``.frq`` files that define allele frequency for all SNPs in the reference. 
    The files may be split per chromosome, similarly to ``bim_file`` argument.
  * [BOTH] ``plink_ld_bin``, required -- path to ``.ld.bin`` files generated by ``bgmg-cli`` as described earlier in this tutorial.
    Note that ``.ld.bin`` files save triplets *(indexA, indexB, r2)*,
    where SNP indices correspond to the reference of SNPs provided to ``bgmg-cli`` during the conversion from plink ``.ld.gz`` files.
	Make sure that ``bim_file`` and ``chr_labels`` arguments are consistent between ``bgmg-cli`` call and ``UGMG_cpp_run_simple`` / ``BGMG_cpp_run_simple`` calls.
  * [BOTH] ``extract``, optional -- file containing SNP rs# to for GWAS SNPs include from the anslysis (fit procedure, QQ plots, and power plots)
  * [BOTH] ``exclude``, optional -- file containing SNP rs# to for GWAS SNPs exclude in the anslysis (fit procedure, QQ plots, and power plots)
  * [BOTH] ``bgmg_shared_library``, required -- path to C/C++ plugin, i.e. ``libbgmg.so`` (linux), ``libbgmg.dylib`` (mac) or ``bgmg.dll`` (windows) 
  * [BOTH] ``bgmg_shared_library_header``, required -- path to ``bgmg_matlab.h`` header

* output data files
  * [BOTH] ``out_file`` -- prefix to output files

* parameters  
  * [BOTH] ``chr_labels`` -- list of chromosome labels to substitute for ``@`` sign in ``bim_file``, ``frq_file``, ``plink_ld_bin`` files.
    It is allowed to set ``chr_labels`` to 1, to run the model only on chromosome 1. 
	However, the same trick will not work for other chromosomes, for the reason described in ``plink_ld_bin``.
  * [UGMG] ``DO_FIT_UGMG``, default ``true`` -- a flag indicating whether MiXeR should perform the fit step; ``false`` is useful to, for example, re-create QQ plots.
  * [BGMG] ``DO_FIT_BGMG``, default ``true`` -- a flag indicating whether MiXeR should perform the fit step; ``false`` is useful to, for example, re-create stratified QQ plots.
  * [BOTH] ``CI_ALPHA``, default ``0.05`` -- significance level for confidence interval estimation; ``nan`` disables uncertainty analysis
  * [BOTH] ``FIT_FULL_MODEL``, default ``true`` -- a flag indicating whether MiXeR should use the full model; setting ``false`` ensures that only fast model is used; this is helpful for diagnostics and debugging
  * [UGMG] ``QQ_PLOT``, default ``true`` -- a flag indicating whether MiXeR should generate QQ plots
  * [UGMG] ``QQ_PLOT_DOWNSCALE``, default ``100`` -- when generating model prediction for QQ plot, MiXeR selects only a small subset of variants; ``QQ_PLOT_DOWNSCALE`` determines which fraction of variants to use. The *data* QQ plot is always based on the entire set of variants available in GWAS, filtered by ``extract`` and ``exclude`` options.
  * [UGMG] ``QQ_PLOT_BINS``, default ``true`` -- a flag indicating whether MiXeR should generate partitioned QQ plots by MAF and LD score
  * [UGMG] ``QQ_PLOT_BINS_DOWNSCALE=``, default ``50`` -- see description in ``QQ_PLOT_DOWNSCALE`` option
  * [UGMG] ``POWER_PLOT``, default ``true`` -- a flag indicating whether MiXeR should generate power curves (proportion of heritability explained as a function of sample size)
  * [UGMG] ``POWER_PLOT_DOWNSCALE``, default ``100`` -- see description in ``QQ_PLOT_DOWNSCALE`` option
  * [BGMG] ``STRATIFIED_QQ_PLOT``, default ``true`` -- a flag indicating whether MiXeR should generate stratified QQ plots (cross-trait enrichment)
  * [BGMG] ``STRATIFIED_QQ_PLOT_DOWNSCALE``, default ``100`` -- see description in ``QQ_PLOT_DOWNSCALE`` option 
  * [BGMG] ``qq_zgrid_lim``, default ``25`` -- defines the range of z-scores for model prediction on stratified QQ plots
  * [BGMG] ``qq_zgrid_step``, default ``0.05`` -- defines the delta step of z-scores for model prediction on stratified QQ plots
  * [BOTH] ``z1max``, default ``nan`` -- enable censoring feature; z-scores exceeding ``z1max`` threshold will contribute to log likelihood via ``cdf(z1max)``.
  * [BGMG] ``z2max``, default ``nan`` -- enable censoring for the second trait;
  * [BOTH] ``randprune_r2``, default ``0.1`` -- threshold for random pruning (procedure that defines weighting of SNPs in log-likelihood cost function and QQ plots)
  * [BOTH] ``randprune_n``, default ``64`` -- number of random pruning iterations
  * [BOTH] ``kmax``, default ``20000`` -- number of samples in full model; it is safe to reduce this number down to ``100`` for all QQ plots, but fit must use large values (recommended ``kmax=20000``)
  * [BOTH] ``max_causal_fraction``, default ``0.03`` -- maximum fraction of causal variants to consider in full model
  * [BOTH] ``cache_tag_r2sum``, default ``false`` -- a flag indicating whether to enable performance optimization that caches z-scores for GWAS SNPs (has no effect on model results)
  * [BOTH] ``r2min``, default ``0.0`` -- lower threshold of ``r2`` in LD structure; all ``r2`` below ``r2min`` will still contribute via infinitesimal model. In addition to ``r2min``, LD structure is truncated by ``plink --r2 --ld-window-r2 0.05`` command.
  * [BOTH] ``SEED``, default ``123`` -- seed for random number generator to use in random pruning and sampling causal variants in *full* model
  * [BOTH] ``THREADS``, default ``-1`` -- how many OMP threads to use for computation (by default uses all available cores) 
  * [BOTH] ``TolX``, default ``1e-2`` -- tolerance of the argument for fminserach stop criteria 
  * [BOTH] ``TolFun``, default ``1e-2`` -- tolerance of the function for fminserach stop criteria

### Memory usage
  
Memory usage of MiXeR largely consists of the following components:

* In-memory storage of LD matrix. This depends on ``plink_ld_bin``, ``extract``/``exclude`` flags and ``r2min`` parameter.
  Typical values are around **10 GB** of memory, for example if one uses all SNPs from LD score reference, and
  ``plink_ld_bin`` files were generated using ``plink r2 gz --ld-window-r2 0.05``.
* Indices of causal variants for sampling in *full* model. This data structure  consumes
  consume ``num_components * kmax * num_snps * max_causal_fraction * 4`` bytes, where ``num_components`` is 1 for univariate and 3 for bivariate,
  ``num_snps`` is the total size of reference (regardless of ``extract`` or ``exclude`` flags).
  With default setting, univariate analysis consumes **22.3 GB** (``kmax=20000``, ``max_causal_fraction=0.03``, ``num_snps=9997231``, ``num_components=1``),
  and bivariate analysis consumes **44.7 GB** (``kmax=20000``, ``max_causal_fraction=0.02``, ``num_snps=9997231``, ``num_components=3``).
* If ``cache_tag_r2sum`` is set to ``true``, MiXeR will consume additional ``num_components * num_gwas_snps * kmax * 4`` bytes,
  where ``num_gwas_snps`` is the number of SNPs with defined z-score, that pass ``extract``/``exclude`` filtering.
  With default settings, this consumes **90.7 GB** per component (``kmax=20000``, ``num_gwas_snps=1217312``),
  therefore this feature is disabled.

## Visualize MiXeR results

Preliminary visualization scripts are available in [vis.py](vis.py) script.
