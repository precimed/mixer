## Contents

* [Introduction](#introduction)
* [Install MiXeR](#install-mixer)
* [Data downloads](#data-downloads)
* [Data preparation](#data-preparation)
* [Run MiXeR](#run-mixer)
* [MiXeR options](#mixer-options)
* [Visualize MiXeR results](#visualize-mixer-results)

## Introduction

Mixer code is generally implemented in Matlab, but some routines were coded in native C or C++ language to give better performance. Therefore, to run MiXeR one needs to either compile C/C++ code, or install per-built binaries which MiXeR depends on. Further details are available in [Install MiXeR](#install-mixer) section.

Input data for MiXeR consists of summary statistics from a GWAS, and a reference panel. MiXeR format for summary statistics is compatible with LD Score Regression (i.e. the ``sumstats.gz`` files), and for those users who are already familiar with ``munge_sumstats.py`` script we recommend to use LD Score Regression pipeline to prepare summary statistics. At the same time, we encourage everyone to take a look [our own pipeline](https://github.com/precimed/python_convert/) for processing summary statistics. For the reference panel we recommend to use 1000 Genomes Phase3 data, pre-processed according to LD Score Regression pipeline, and available for download from LDSC website. Further details are given in [Data downloads](#data-downloads) and [Data preparation](#data-preparation) sections.

Once you have all input data in MiXeR-compatible format you may proceed with running univariate analysis ([UGMG_cpp_run_simple.m](UGMG_cpp_run_simple.m) script) and cross-trait analysis ([BGMG_cpp_run_simple.m](BGMG_cpp_run_simple.m) script). The results will be saved as ``.json`` files. To visualize the results we provide a script in python, but we encourage users to write their own scripts that process the results. Further details are given in [Run MiXeR](#run-mixer), [MiXeR options](#mixer-options) and [Visualize MiXeR results](#visualize-mixer-results) sections.

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
  * change active folder to ``MIXER_ROOT``
  * add ``$MIXER_ROOT/lib`` to ``LD_LIBRARY_PATH`` (so that matlab can find boost libraries from MiXeR distribution package):
    ```
    export "LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH"
    ```
  * start matlab and execute ``test_mixer_plugin`` command:
    ```
    matlab -nodisplay -nosplash -nodesktop -r "test_mixer_plugin; exit;"
    ```
    
  * check that ``test_mixer_plugin`` have created a new file named ``test_mixer_plugin.bgmglib.log`` containing the following lines:
    ```
    20181208 20:16:56.281717	============= new session =============
    20181208 20:16:56.281717	=mixer plugin setup success
    ```
    
This procedure was tested on
  * ``CentOS release 6.10`` (on UCSD Comet cluster):
    * ``Matlab R2018a (9.4.0.813654) 64-bit (glnxa64)`` - all OK
  * ``CentOS release 6.9`` (on UiO Abel cluster):
    * ``Matlab R2018a (9.4.0.813654) 64-bit (glnxa64)`` - all OK
    * ``Matlab R2017a (9.2.0.556344) 64-bit (glnxa64)`` - all works but matlab crashes upon exit
    * ``Matlab R2016a (9.0.0.341360) 64-bit (glnxa64)`` - all works but matlab crashes upon exit
  * ``Debian GNU/Linux 9`` (on LISA supercomputer at SURFsara)
    * ``Matlab R2017b (9.3.0.713579) 64-bit (glnxa64)`` - all works but matlab crashes upon exit
    * ``Matlab R2016b (9.1.0.441655) 64-bit (glnxa64)`` - all works but matlab crashes upon exit

### Install on Windows using pre-built binaries

* Download "Windows_x64.7z" file from the latest MiXeR release (https://github.com/precimed/mixer/releases)
* Extract "Windows_x64.7z" to a new folder. Below we refer to this folder as ``MIXER_ROOT``.
* Test that MiXeR executable runs smoothly, as described in [Install on Linux using pre-built binaries](#install-on-linux-using-pre-built-binaries) section
* Test that MiXeR C++ plugin is loaded correctly, as described in [Install on Linux using pre-built binaries](#install-on-linux-using-pre-built-binaries) section

### Build from source - Linux

We build release package on Abel supercomputer.

1. Download [boost/1_49_0](https://www.boost.org/doc/libs/1_49_0/more/getting_started/unix-variants.html), and compile it
   ```
   module load gcc/4.7.2
   ./bootstrap.sh --with-libraries=program_options,filesystem,system,date_time && ./b2 --j12
   ```
2. Clone MiXeR repository, and compile it
   ```
   git clone --recurse-submodules -j8 git://github.com/precimed/mixer.git && cd mixer/src
   module purge && module load cmake/3.7.1 gcc/4.7.2  #  cmake/3.7.1 gcc/4.7.2
   mkdir build && cd build && cmake .. -DBOOST_ROOT=/usit/abel/u1/oleksanf/boost_1_49_0  
   make -j12 && cd ..
   ```
3. Compile ``bgmg-cli`` separately (usenother combination of modules)
   ```
   module purge
   module load cmake/3.7.1 gcc/4.9.0 openmpi.gnu/1.8.1 icu/49.1.2
   mkdir build2 && cd build2 && cmake .. -DBOOST_ROOT=/usit/abel/u1/oleksanf/boost_1_49_0
   make -j12 && cd ..
   ```

The specific combination of ``gcc`` and ``boost`` versions is important - it turned out to be compatible with matlab.
Further notes are available in [src/README.md](src/README.md).
		
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
  tar -xzvf 1000G_Phase3_plinkfiles.tgz
  bzip2 -d w_hm3.snplist.bz2
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
       --out LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p05_SNPwind50k
    ```
    This step takes about 1 hour (here and below all times are measured on a server with 24 physical cores).
  * Run ``bgmg-cli`` to convert plink output into a binary format. The following command must be run once for each chromosome. 
    Note that ``--bim`` argument must be the same, i.e. point to ``1000G.EUR.QC.@.bim`` regardless of the actual chromosome that you use in ``--plink-ld`` and ``--out``.
    ```
    bin/bgmg-cli \
       --bim LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
       --plink-ld LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p05_SNPwind50k.ld.gz \
       --out LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p05_SNPwind50k.ld.bin
    ```
    The conversion is single-threaded, and takes about 10 minutes for the longest chromosome.
    The output is written to ``LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p05_SNPwind50k.ld.bin`` file,
    and log details into ``LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr_label>.p05_SNPwind50k.ld.bin.bgmglib.log`` file:
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
    ```
  * Save the list of dbSNP rs# into a separate file called ``w_hm3.justrs``:
    ```
    cut -f1 w_hm3.snplist | tail -n +2 > w_hm3.justrs
    ```

## Run MiXeR

### Univariate analysis

To start univariate analysis one should start matlab, configure parameters as listed below, and execute [UGMG_cpp_run_simple](UGMG_cpp_run_simple.m) script.
The recommended way is to call ``cd $MIXER_ROOT && matlab -nodisplay -nosplash -nodesktop -r "UGMG_params_fit; UGMG_cpp_run_simple; exit;"``, where ``UGMG_params_fit.m`` is a file that looks as follows:
```
trait1_file='SSGAC_EDU_2018_no23andMe_noMHC.sumstats.gz';
out_file='SSGAC_EDU_2018_no23andMe_noMHC.fit';
bim_file='LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim';
frq_file='LDSR/1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.@.frq';
plink_ld_bin='LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin'; chr_labels = 1:22;
init_from_params_file=''; extract='LDSR/w_hm3.justrs';
bgmg_shared_library='lib/libbgmg.so';
bgmg_shared_library_header='bgmg_matlab.h';
kmax=20000; max_causal_fraction=0.03; z1max=nan; z2max=nan;
cache_tag_r2sum=0; r2min=0.0;
randprune_r2=0.1; randprune_n=64;
DO_FIT_UGMG=1; FIT_FULL_MODEL=1; CI_ALPHA=0.05; SEED=123;
POWER_PLOT=0; POWER_PLOT_DOWNSCALE=100;
QQ_PLOT=0; QQ_PLOT_DOWNSCALE=100;
QQ_PLOT_BINS=0; QQ_PLOT_BINS_DOWNSCALE=50;
```
The results will be saved ``<out_file>.json`` file.

The above parameters imply that the model is fitted on HapMap3 SNPs, as specified by ``extract='w_hm3.justrs'``.
To produce QQ plots on the entire set of SNPs one may change parameters as follows:
```
DO_FIT_UGMG=0; kmax=100; extract=''; 
POWER_PLOT=1; QQ_PLOT=1; QQ_PLOT_BINS=1;
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
init_from_params_file=''; extract='LDSR/w_hm3.justrs';
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

To produce QQ plots on the entire set of SNPs one may change parameters as follows:
```
DO_FIT_BGMG=0; kmax=100; extract='';
STRATIFIED_QQ_PLOT=1;
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

## MiXeR options

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

### MiXeR results format

All MiXeR results are stored in a single ``.json`` file.

Results of univariate analysis:
* ``['univariate'][0]['params'][<parameter>]`` - point estimates of model parameters. Here ``<parameter>`` can be one of the following: ``pi_vec`` (polygenicity), ``sig2_beta`` (variance of causal effect sizes), ``sig2_zero`` (variance distortion)
* ``['univariate'][0]['ci'][<parameter>][<measure>]`` - uncertainty of parameter estimates. Here ``<parameter>`` can be one of the following: ``h2`` (heritability), ``pi_vec``, ``sig2_beta``, ``sig2_zero`` -  as described above; ``<measure>`` can be one of the following: ``point_estimate`` (point estimate),  ``se`` (standard error), ``lower`` and ``upper`` (lower and upper bound of confidence interval, at significance level defined by ``CI_ALPHA`` option);
* ``['univariate'][0]['power_plot_data']['power_nvec']`` and ``['univariate'][0]['power_plot_data']['power_svec']`` - power plot; ``power_svec`` gives fraction of heritability explain for a given sample size (``power_nvec``).
* ``['univariate'][0]['qq_plot_data'][<variable>]`` - data and model QQ plots. ``<variable>`` can be one of the following: 
``hv_logp`` - observed log(p-values), ``data_logpvec`` -- expected log(p-values) for the data; ``model_logpvec`` -- expected log(p-values) for the model;
* ``['univariate'][0]['qq_plot_bins_data'][<index>]`` - a 3x3 matrix of QQ plots, partitioned by MAF and LD score; ``<index>`` can take values 0 to 8; each QQ plot follows the format defined above; the ranges of MAF and LD score for each bin is specified in ``['univariate'][0]['qq_plot_bins_data'][<index>]['title']``.

Results of bivariate analysis:
* ``['bivariate']['params'][<parameter>]`` - point estimates of model parameters. Here ``<parameter>`` can be one of the following: ``pi_vec`` - polygenicty of each component (fraction of variants specific to the first trait, specific to the second trait, and shared across traits); ``rho_beta`` - correlation of effect sizes within each component (first two values are zeros); ``rho_zero`` - covariance inflation parameter; ``sig2_beta`` - 2x3 matrix, variance of effect sizes for each trait and within each component; ``sig2_zero`` - variance distortion in each trait. 
* ``['bivariate']['ci'][<parameter>][<measure>]`` - uncertainty of parameter estimates. Here ``<parameter>`` can be one of the following: ``h2_T1``, ``h2_T2`` - heritability of the first and the second traits; ``pi1u``, ``pi2u`` total polygenicity of the first and of the second trait; ``pi_vec_C1``, ``pi_vec_C2``, ``pi_vec_C3`` - polygenicity of the three components in the model; ``rg`` - genetic correlation; ``rho_beta`` - correlation of effect sizes within shared polygenic component; ``rho_zero`` - covariance inflation parameter; ``sig2_beta_T1``, ``sig2_beta_T2`` - variance of effect sizes in the first and in the second trait; ``sig2_zero_T1``, ``sig2_zero_T1`` - variance distortion in the first and in the second trait.)
* ``['bivariate']['stratified_qq_plot_fit_data'][<trait>][<index>]`` - data for stratified QQ plots. Here ``<trait>`` can be either ``'trait1'`` or ``'trait2'``; ``<index>`` can be ``0``, ``1``, ``2`` or ``3``. Stratified QQ plots are made both ways (i.e. ``trait1|trait2``, and ``trait2||trait1``; hense the ``<trait>`` defines which trait is primary (so that stratified QQ plots are conditioned on the other trait. ``<index>`` equal to ``0`` means a stratum of all SNPs, ``1``, ``2`` and ``3`` are increased levels of association on the secondary trait. Each QQ plot has format defined above (see results of univariate analysis).

