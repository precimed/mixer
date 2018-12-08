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

* Reference pannel
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

## Visualize MiXeR results

Preliminary visualization scripts are available in [vis.py](vis.py) script.
