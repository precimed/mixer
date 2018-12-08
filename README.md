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
The results will be saved to files named according to ``out_file`` parameter.

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
The results will be saved to files named according to ``out_file`` parameter.

To test fitted  parameters on the entire set of SNPs one may change parameters as follows:
```
DO_FIT_BGMG=0; kmax=100; extract='';
init_from_params_file='PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.fit.params.mat';
out_file='PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.test';
```

## Visualize MiXeR results

# (legacy stuff) UGMG and BGMG - Univariate and Bivariate Gaussian Mixtures for GWAS summary statistics

`UGMG` model allows you to estimate polygenicity (e.i. proportion of causal variants) of a complex trait given GWAS summary statistics. 
`BGMG` model allows you to quantify polygenic overlap (e.i. proportion of causal variants shared between two complex traits) given GWAS summary statistics. As of today both methods published on biorxiv:
* D. Holland et al., [Estimating phenotypic polygenicity and causal effect size variance from GWAS summary statistics while accounting for inflation due to cryptic relatedness](https://www.biorxiv.org/content/early/2017/06/23/133132)
* O. Frei et al., [Bivariate Gaussian Mixture Model of GWAS (BGMG) quantifies polygenic overlap between complex traits beyond genetic correlation](https://www.biorxiv.org/content/early/2017/12/27/240275)

## Requirements

1. ``Matlab`` (tested with ``8.6.0.267246 (R2015b)``, may work with other versions as well)
2. ``python 3 > version >= 2.7`` --- to convert summary stats into ``BGMG``-compatible format. 

## Getting Started

To download `BGMG` you should clone these two repositories:
```  
git clone https://github.com/precimed/bgmg.git            # or git@github.com:precimed/BGMG.git
git clone https://github.com/precimed/python_convert.git  # or git@github.com:precimed/python_convert.git
```

``BGMG`` repository contains the actual matlab code of the models (both UGMG and BGMG).
``python_convert`` repository contains supplementary code that loads summary statistics into ``BGMG``-compatible format.

You also need to download reference data:
```
wget "https://www.dropbox.com/s/5cvqzaayg1tn5u0/all_chromosomes_multiallelic_replaced_with_AT.ref.gz?dl=1" -O all_chromosomes_multiallelic_replaced_with_AT.ref.gz
wget "https://www.dropbox.com/s/u3nhdznhcun7tok/1kG_phase3_EUR_11015883_reference_holland.mat?dl=1" -O 1kG_phase3_EUR_11015883_reference_holland.mat
wget "https://www.dropbox.com/s/nrhh96w8tzavcus/ldmat_p8_BPwind10M_n503.mat?dl=1" -O ldmat_p8_BPwind10M_n503.mat
wget "https://www.dropbox.com/s/j08f848raabcmcf/all_chromosomes.ref.gz?dl=1" -O all_chromosomes.ref.gz
wget "https://www.dropbox.com/s/tcd4hmu71ks0xny/defvec_1kG_phase3_EUR.mat?dl=0" -O defvec_1kG_phase3_EUR.mat
```
``1kG_phase3_EUR_11015883_reference_holland.mat`` and ``ldmat_p8_BPwind10M_n503.mat`` contains information about linkage disequilibrium structure, extracted for EUR population from 1000 genomes project, phase 3. Additional details are given further below.
``all_chromosomes_multiallelic_replaced_with_AT.ref.gz`` contains a list of 11015833 variants that we use to fit UGMG and BGMG model, with additional hack to replace multiallelic SNPs with ``A/T`` to workaround limitations of ``sumstats.py`` script. The original reference file is also enclosed (all_chromosomes.ref.gz). It was produced after basic QC steps, applied to all SNPs from 1000 genomes project, phase 3. 

For further information about the reference file please refer to [out wiki pages](https://github.com/precimed/BGMG/wiki/How-to-create-a-reference-file).
For additional data files, including simulation data, look [here](https://www.dropbox.com/sh/kpwj02vxlbm2fh9/AABsep2F2xJGRK1iX198L5Bja?dl=0).

## Prepare summary statistics

First step is to convert summary statistics into `BGMG`-compatible format, using ``sumstats.py`` tool from ``python_convert`` repository. The following example shows how to convert Schizophrenia ([PGC, 2014, scz2.snp.results.txt.gz](https://www.med.unc.edu/pgc/results-and-downloads)) and Years of Education ([SSGAC, 2016, EduYears_Main.txt](https://www.thessgac.org/data)).

```
python sumstats.py csv \
	--sumstats scz2.snp.results.txt.gz \
	--ncase-val 35476 \
	--ncontrol-val 46839 \
	--auto  \
	--chr hg19chrc \
	--out PGC_SCZ_2014.csv && gzip PGC_SCZ_2014.csv \

python sumstats.py csv \
	--sumstats EduYears_Main.txt.gz \
	--auto  \
	--n-val 293723 \
	--out SSGAC_EDU_2016.csv && gzip SSGAC_EDU_2016.csv \
  
python sumstats.py mat --sumstats PGC_SCZ_2014.csv.gz   --ref all_chromosomes_multiallelic_replaced_with_AT.ref.gz --out PGC_SCZ_2014.mat
python sumstats.py mat --sumstats SSGAC_EDU_2016.csv.gz --ref all_chromosomes_multiallelic_replaced_with_AT.ref.gz --out SSGAC_EDU_2016.mat
```
For more information about ``sumstats.py`` utility see ``python_convert`` repository,
or type ``sumstats.py --help`` or ``sumstats.py <command> --help`` where ``<command>`` is either ``csv`` or ``mat``.

This results in a separate ``.mat`` file for each trait.
Each file contains two variables, ``zvec`` and ``nvec`` ---
a signed test statistic (`z-score`) and effective sample size (`nvec`) for each variant. 
``zvec`` and ``nvec`` will be a column vector of length `N` (number of variants from the reference file, ``all_chromosomes.ref.gz``).
The order of the variants will match the reference file.
``zvec`` or ``nvec`` might contain undefined values, ``nan``.
In this case the analysis will be limited to a subset of variants where ``zvec`` and ``nvec`` are defined.

## Running univariate model (UGMG)

To run UGMG model you need to ``cd`` to the root of BGMG repository. Then you may run UGMG from console as follows.
```
matlab -nodisplay -nosplash -nodesktop -r "trait1_file='PGC_SCZ_2014.mat'; reference_file='1kG_phase3_EUR_11015883_reference_holland.mat'; QQ_PLOT_FIT=true; LDmat_file='ldmat_p8_BPwind10M_n503.mat'; LDmat_file_variable='LDmat'; out_file='PGC_SCZ_2014.result'; defvec_files={'defvec_1kG_phase3_EUR.mat'}; BGMG_run; exit;"
```
Alternatively, you may open matlab, create variables ``trait1_file``, ``reference_file``, ``out_file`` with values as in the script above, and then execute ``BGMG_run.m`` script. The meaning of the parameters is as follows:
* ``trait1_file`` points to the output of ``sumstats.py``, e.i. summary stats converted to matlab format. It must be a filename, without path, of the ``.mat`` file containing ``zvec`` of the trait to be analyzed.
* ``reference_data`` should be set to absolute or relative path of ``B1kG_phase3_EUR_11015883_reference_holland.mat``
* ``defvec_files`` specifies a cell array of files with defvec. Each defvec is a binary column-vector of length ``11015883`` (e.i. number of variants in the reference) with ``0`` for variants to exclude from the analysis (for example MHC, low maf, etc).
* ``out_file`` specifies the name of the resulting file

## Interpret univariate results (UGMG)

After execution you should see approximately the following log file,
and a ``.pdf`` file with the following [QQ plot](https://github.com/precimed/BGMG/blob/master/PGC_SCZ_2014.result.pdf).

```
>> BGMG_run
trait1: PGC_SCZ_2014.mat
Loading reference data from BGMG_reference_data_11Msnps_2018_02_12.mat...
11015833 SNPs in the template
6610991 SNPs left after filtering missing zvec and nvec values, trait PGC_SCZ_2014.mat
6610991 SNPs left after filtering missing values (maf, tld, etc)
6506468 SNPs left after filtering large LD blocks (<= 2000.00)
6476539 SNPs left after filtering high LD SNPs (<= 600.00)
6464140 SNPs left after filtering MHC
6269232 SNPs left after filtering low MAF SNPs (>= 0.010)
Perform 10 iterations of random pruning ..........done.
Effective number of SNPs on each iteration of random pruning:
ans =
     1484784     1484826     1484320     1484902     1484690     1484855     1484416     1484786     1484682     1484857
Trait1  : fit infinitesimal model to find initial sig2_zero
Univariate: pi_vec=1.000e+00, sig2_beta^2=1.238e-05, sig2_zero^2=1.607, h2=29.549, cost=3.064280e+06, nsnp=3678219
Univariate: pi_vec=1.000e+00, sig2_beta^2=1.238e-05, sig2_zero^2=1.646, h2=29.549, cost=3.067946e+06, nsnp=3678219
...
Trait1  : fit pi_vec and sig2_beta, constrained on sig2_zero and h2
Univariate: pi_vec=1.000e-02, sig2_beta^2=1.780e-05, sig2_zero^2=1.172, h2=0.425, cost=2.324409e+06, nsnp=3678219
Univariate: pi_vec=7.964e-03, sig2_beta^2=2.235e-05, sig2_zero^2=1.172, h2=0.425, cost=2.324302e+06, nsnp=3678219
...
Trait1  : final unconstrained optimization
Univariate: pi_vec=3.168e-03, sig2_beta^2=5.619e-05, sig2_zero^2=1.172, h2=0.425, cost=2.323974e+06, nsnp=3678219
Univariate: pi_vec=2.378e-03, sig2_beta^2=5.619e-05, sig2_zero^2=1.172, h2=0.319, cost=2.324707e+06, nsnp=3678219
...
Trait1  : optimization with poisson cost function
Univariate: pi_vec=3.122e-03, sig2_beta^2=5.746e-05, sig2_zero^2=1.171, h2=0.428, cost=2.323818e+06, nsnp=3678219
Univariate: pi_vec=3.122e-03, sig2_beta^2=5.746e-05, sig2_zero^2=1.171, h2=0.428, cost=2.323818e+06, nsnp=3678219
...
done.
Estimate cumulated distribution function for Z scores...
Finish 313461 SNPs out of 6269232 - Poisson grid 0:0.09:2.16 (25 nodes)
Finish 626922 SNPs out of 6269232 - Poisson grid 0:0.10:2.45 (25 nodes)
...
Univariate: pi_vec=3.039e-03, sig2_beta^2=5.989e-05, sig2_zero^2=1.164, h2=0.434, cost=1.026881e+07, nsnp=6269232
Result saved to PGC_SCZ_2014.result.pdf
```

## Running bivariate model (BGMG)

TBD.

## Interpret bivariate results (BGMG)

TBD.

## Authors

* Dominic Holland (Center for Multimodal Imaging and Genetics, University of California at San Diego)
* Oleksandr Frei (NORMENT, University of Oslo)
* Anders M. Dale (Center for Multimodal Imaging and Genetics, University of California at San Diego)

## Legacy instructinos on how to create a reference

Format of the reference file
----------------------------

***
TBD: include hvec variable in the reference. Discuss implications of GCTA maf model (LDAK alpha=-1) vs BGMG maf model (LDAK alpha=0).
***

[This](https://www.dropbox.com/s/7trec3ma6u1m6uy/HAPGEN_EUR_100K_11015883_reference_holland.mat?dl=0) 
is an example of a reference file, compatible with UGMG and BGMG model.
Reference file is based on certain genotype reference (for example from 1000 genomes project).

Reference file contains the following variables.
* `chrnumvec` - chromosome labels, integer column vector, one value for each variant in the reference, coded from 1 to 22
* `posvec` - base pair positions, integer column vector, one value for each variant in the reference, unity-based (as ENSEMBL)
* `mafvec` - minor allele frequencies, double-precision column vector, one value for each variant in the reference.
* `total_het` - double-precision scalar value. The value depends on MAF model; for GCTA model (LDAK alpha = -1) total_het must be equal to the number of genotyped SNPs. For BGMG model (LDAK alpha = 0) total heterozigosity must be calculated as a sum of `2*maf*(1-maf)` across all genotyped variants.
* `hvec` - double-precision column vector, one value for each variant in the reference. The value depends on MAF model; for GCTA model (LDAK alpha = -1) all hvec values must must equal to 1. For BGMG model (LDAK alpha = 0) ``hvec=2*mafvec.*(2-mafvec)``.
* `w_ld` - LD score, double-precision column vector, one value for each variant in the reference. LD scores must be calculated towards variants within the reference.
* `sum_r2` - LD scores, double-precision matrix, rows correspond to variants in the reference, columns correspond to bins of allelic correlation (for example 0.05 to 0.25, 0.25 to 0.50, 0.50 to 0.75, and 0.75 to 1.00). 
* `sum_r2_biased` - same format as `sum_r2`, but based on allelic correlation that was not corrected for bias
* `sum_r4_biased` - same format as `sum_r2`, but based on 4-th power of allelic correlation that was not corrected for bias

Note that `sum_r2`, `sum_r2_biased`, `sum_r4_biased` depend on MAF model. For BGMG model (LDAK alpha = 0) one should use ``--per-allele`` flag in ``ldsc.py --ld``; otherwise, for GCTA model (LDAK alpha = -1), there calculation must run without ``--per-allele``.

`posvec` and `chrnumvec` are optional, but we recommend to include them to keep track of what variants are used in the analysis.
We don't recommend to save list of variant IDs, such as RS numbers, because of performance limitations in matlab (saving long cell arrays is too slow).

Tools required to gather reference file
---------------------------------------

1. Modified version of LDSC software:
```
git clone https://github.com/ofrei/ldsc.git
```
2. `sumstats.py` script:
```
git clone https://github.com/precimed/python_convert.git
```
3. [plink](https://www.cog-genomics.org/plink2), to calculate allele frequencies and, optionally, gather pairwise LD correlations
4. You also need reference genotypes, split per chromosome. We recommend to download [these](https://www.dropbox.com/s/sgzhmc1q1me6zji/1kG_phase3_EUR_11015883.tar.gz?dl=0) data, derived from 1kG phase 3.

Scripts
-------
Change ``LDSC_PY``, ``SUMSTATS_PY``, ``BFILE``, ``REF``, ``OUT_FOLDER`` to match your path.
Optionally, set ``PYTHON`` and ``PLINK`` to the full path to your python and plink executables.
Step *calculate LD structure* may take considerably long time,
and if you have access to an HPC cluster we recommend to submit commands from ``commands.txt`` and them all of them in parallel.
```
# [BASH] Setup code and data location 
SUMSTATS_PY=/mnt/h/GitHub/python_convert/sumstats.py
LDSC_PY=/mnt/h/GitHub/ldsc/ldsc.py
PYTHON=python
PLINK=plink
BFILE="/mnt/h/NORSTORE/oleksanf/1kG_phase3_EUR_11015883"
REF=/mnt/h/Dropbox/shared/BGMG/all_chromosomes.ref.gz
OUT_FOLDER="/mnt/h/NORSTORE/oleksanf/1kG_phase3_EUR_11015883_reference"


R2_BINS_MIN="0.05 0.25 0.50 0.75"
R2_BINS_MAX="0.25 0.50 0.75 1.00"
LD_WIND_SNPS=50000

# [PYTHON] Create ref file (e.i. combined bim file with header) from chromosome-split genotypes using python script:
import pandas as pd
import os
df = pd.concat([pd.read_table('/mnt/h/NORSTORE/oleksanf/1kG_phase3_EUR_11015883/chr{}.bim'.format(chri), header=None, names=['CHR', 'SNP', 'GP', 'BP', 'A1', 'A2']) for chri in range(1, 23)]).to_csv('all_chromosomes.ref.gz', sep='\t', compression='gzip', index=False)

# [BASH] calculate allele frequencies and convert them to matlab
for CHR in {1..22}; do $PLINK --bfile $BFILE/chr$CHR --freq --out $OUT_FOLDER/chr$CHR; done
$PYTHON $SUMSTATS_PY frq-to-mat --ref $REF --force --frq $OUT_FOLDER/chr@.frq --out $OUT_FOLDER/mafvec.mat

# [BASH] load reference file into matlab
$PYTHON $SUMSTATS_PY ref-to-mat --ref $REF --out $OUT_FOLDER/all_chromosomes.mat --force --numeric-only

# [BASH] calculate LD structure
rm commands.txt
for CHR in {1..22} ; do echo "$PYTHON $LDSC_PY --l2           --ld-wind-snps $LD_WIND_SNPS --bfile $BFILE/chr$CHR                                              --out $OUT_FOLDER/w_ld_chr$CHR" >> commands.txt; done
for CHR in {1..22} ; do echo "$PYTHON $LDSC_PY --l2           --ld-wind-snps $LD_WIND_SNPS --bfile $BFILE/chr$CHR --r2-min $R2_BINS_MIN --r2-max $R2_BINS_MAX  --out $OUT_FOLDER/ref_ld_r2_unbias_chr"$CHR >> commands.txt; done
for CHR in {1..22} ; do echo "$PYTHON $LDSC_PY --l2 --bias-r2 --ld-wind-snps $LD_WIND_SNPS --bfile $BFILE/chr$CHR --r2-min $R2_BINS_MIN --r2-max $R2_BINS_MAX  --out $OUT_FOLDER/ref_ld_r2_biased_chr"$CHR >> commands.txt; done
for CHR in {1..22} ; do echo "$PYTHON $LDSC_PY --l4 --bias-r2 --ld-wind-snps $LD_WIND_SNPS --bfile $BFILE/chr$CHR --r2-min $R2_BINS_MIN --r2-max $R2_BINS_MAX  --out $OUT_FOLDER/ref_ld_r4_biased_chr"$CHR >> commands.txt; done
cat commands.txt | parallel -j12

$PYTHON $SUMSTATS_PY ldsc-to-mat --ref $REF --force --ldscore $OUT_FOLDER/w_ld_chr@.l2.ldscore.gz --out $OUT_FOLDER/w_ld.l2.mat
for R2BIN in {1..4} ; do $PYTHON $SUMSTATS_PY ldsc-to-mat --ref $REF --force --ldscore $OUT_FOLDER/ref_ld_r2_unbias_chr@-r2bin-$R2BIN.l2.ldscore.gz --out $OUT_FOLDER/ref_ld_r2_unbias-r2bin-$R2BIN.l2.mat ; done
for R2BIN in {1..4} ; do $PYTHON $SUMSTATS_PY ldsc-to-mat --ref $REF --force --ldscore $OUT_FOLDER/ref_ld_r2_biased_chr@-r2bin-$R2BIN.l2.ldscore.gz --out $OUT_FOLDER/ref_ld_r2_biased-r2bin-$R2BIN.l2.mat ; done
for R2BIN in {1..4} ; do $PYTHON $SUMSTATS_PY ldsc-to-mat --ref $REF --force --ldscore $OUT_FOLDER/ref_ld_r4_biased_chr@-r2bin-$R2BIN.l4.ldscore.gz --out $OUT_FOLDER/ref_ld_r4_biased-r2bin-$R2BIN.l4.mat ; done

# [MATLAB] Combine all files together into a single mat file.
cd $OUT_FOLDER
load all_chromosomes.mat; posvec = BP; chrnumvec=CHR;
load mafvec.mat; total_het=sum(2*mafvec.*(1-mafvec));
load w_ld.l2.mat; w_ld = annomat;
num_r2bins=4;sum_r2=nan(length(mafvec), num_r2bins); sum_r2_biased = sum_r2; sum_r4_biased = sum_r2;
for r2bin=1:num_r2bins
data=load(sprintf('ref_ld_r2_unbias-r2bin-%i.l2.mat', r2bin)); sum_r2(:, r2bin) = data.annomat;
data=load(sprintf('ref_ld_r2_biased-r2bin-%i.l2.mat', r2bin)); sum_r2_biased(:, r2bin) = data.annomat;
data=load(sprintf('ref_ld_r4_biased-r2bin-%i.l4.mat', r2bin)); sum_r4_biased(:, r2bin) = data.annomat;
end
save('reference.mat', 'posvec', 'chrnumvec', 'total_het', 'mafvec', 'w_ld', 'sum_r2', 'sum_r2_biased', 'sum_r4_biased');
```

Remark about genotyped variants, and variants in the reference 
--------------------------------------------------------------
Generally, reference file might represent a subset of the genotyped variants (for example HapMap3 variants as opposite to 1kG phase3 genotypes). You reduce you reference to a subset of genotyped variants. In this case is important to keep in mind that:
* ``w_ld`` must be calculated within the reference, e.i. ignore genotype variants not in the reference
* All ``ref_ld`` (biased_r2, biased_r4, unbias_r2) must be calculated towards all genotyped variants
* ``total_het`` must be calculated across all genotyped variants

[Optional] Binary pairwise LD matrix for random pruning
-------------------------------------------------------

This step is optional, and provides an alternative to weighting SNPs by inverse ``w_ld`` score.
The alternative is based on *random pruning* technique, 
which randomly selects a subsets of pseudo-independent variants
(*pseudo* means that their pairwise LD r2 will be below certain threshold).

To calculate pairwise LD matrix:
```
[BASH] Use plink to calculate pairwise LD
seq 1 22 | parallel -j12 "$PLINK --memory 10240 --threads 2 --bfile $BFILE/chr{} --r2 gz --ld-window-kb 10000 --ld-window 999999 --ld-window-r2 0.8 --out $OUT_FOLDER/ldmat_p8_BPwind10M_n503_chr{}"

[BASH] Load plink results into matlab
MAKE_LD_MATRIX_PY=/mnt/h/GitHub/python_convert/make_ld_matrix/make_ld_matrix.py
seq 1 22 | parallel -j12 "$PYTHON $MAKE_LD_MATRIX_PY --ldfile $OUT_FOLDER/ldmat_p8_BPwind10M_n503_chr{}.ld.gz --ref $REF --savemat $OUT_FOLDER/ldmat_p8_BPwind10M_n503_chr{}.mat"

[MATLAB] Merge files together
cd $OUT_FOLDER
outname = 'ldmat_p8_BPwind10M_n503.mat';
for chri=1:22, fprintf('%i…\n', chri);
    df = load(sprintf('ldmat_p8_BPwind10M_n503_chr%i.mat', chri), 'id1', 'id2', 'nsnp');
    if chri==1, id1=df.id1; id2=df.id2; nsnp=df.nsnp; else id1 = [id1; df.id1]; id2 = [id2; df.id2]; end;
end;
LDmat = sparse(double(id1),double(id2),true,double(nsnp),double(nsnp));
LDmat = LDmat | speye(double(nsnp));LDmat = LDmat | (LDmat - LDmat');
save(outname, 'LDmat', '-v7.3');clear;
```
The resulting matrix can be passed to ``BGMG_run.m`` via ``LDmat_file='ldmat_p8_BPwind10M_n503.mat'; LDmat_file_variable='LDmat';`` flags. In this case ``w_ld`` data from the reference file will not be used.



