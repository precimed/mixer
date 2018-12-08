Introduction
============

Mixer code is generally implemented in Matlab, but some routines were coded in native C or C++ language to give better performance. Therefore, to run MiXeR one needs to either compile C/C++ code, or install per-built binaries which MiXeR depends on. Further details are available in "Install MiXeR" section.

Input data for MiXeR consists of summary statistics from a GWAS, and a reference panel. We generally recommend to use LD Score Regression pipeline to prepare summary statistics. For the reference panel we recommend to use 1000 Genomes Phase3 data, pre-processed according to LD Score Regression pipeline, and available for download from LDSC website. Further details are given in "Data downloads" and "Data preparation" sections.

Once you have all input data in MiXeR-compatible format you may proceed with running univariate analysis (``UGMG_cpp_run_simple.m`` script) and cross-trait analysis (``BGMG_cpp_run_simple.m`` script). The results will be saved as ``.json`` files. To visualize the results we provide a script in python, but we encourage users to write their own scripts that process the results. Further details are given in "Run MiXeR" and "Visualize MiXeR results" sections.

If you encounter an issue, or have further questions, please create a new issue ticket on github.com/precimed/mixer .

If you use MiXeR software for your research publication, please cite the following paper(s):

* for univariate analysis: D. Holland et al., [Estimating phenotypic polygenicity and causal effect size variance from GWAS summary statistics while accounting for inflation due to cryptic relatedness](https://www.biorxiv.org/content/early/2017/06/23/133132) 
* for cross-trait analysis: O.Frei et al., Bivariate causal mixture model quantifies polygenic overlap between complex traits beyond genetic correlation, bioXriv, doi: https://doi.org/10.1101/240275 

The MiXeR software may not be used for commercial purpose or in medical applications.
We encourage all users to familiarize themselves with US patent https://www.google.no/patents/US20150356243 "Systems and methods for identifying polymorphisms".

Install MiXeR
=============

Prerequisites
-------------

* MiXeR was tested on Linux (TBD: abel version, lisa version) and Windows 10 operating systems (TBD: build)
* Matlab (TBD version) (for Windows 10); Matlab (TBD version) (for Linux).
  Other versions may work as well, but are not guarantied to be compatible with pre-built MiXeR binaries (C/C++ code).
* Python 2.7 (for LD score regression)
* Python 3.5 (for MiXeR results visualization)

Hardware requirements
---------------------

MiXeR software is very CPU and memory intensive. 
Minimal memory requirement is to have 61.5 GB of RAM available to MiXeR.
MiXeR efficiently uses multiple CPUs. We recommend to run MiXeR on a system with 16 physical cores.
When use MiXeR on a cluster, we recommend to assign the whole node to each MiXeR run.

Install on Linux (use pre-built binaries)
-----------------------------------------

* Download MiXeR software from github (either 'git clone' from command line, or using the 'download' button on the webpage). You need .m scrits in the root of the repository, and ``DERIVESTsuite`` folder. We refer to this folder as ``MIXER_ROOT``.
* Download pre-built binaries for Linux (see github releases)
* Open command line and  ``cd MIXER_ROOT``. Test that MiXeR C++ plugin is loaded correctly.
* Test that MiXeR binary runs smoothly.
	
Install on Windows (use pre-built binaries)
-------------------------------------------

* Download MiXeR software from github (either 'git clone' from command line, or using the 'download' button on the webpage). You need .m scrits in the root of the repository, and ``DERIVESTsuite`` folder.
* Download pre-built binaries for Windows (see github releases)

Build from source - Linux
-------------------------

TBD. Preliminary notes are available in src/README.md.
		
Build from source - Windows
---------------------------

TBD. Preliminary notes are available in src/README.md.


Data downloads
==============

- LDSR data
- Summary statistics

Data preparation
================

	Step1. Generate LD r2 correlations using plink
	Step2. Convert LD r2 correlations into binary MiXeR format

Run MiXeR
=========

Visualize MiXeR results
=======================


# UGMG and BGMG - Univariate and Bivariate Gaussian Mixtures for GWAS summary statistics

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
