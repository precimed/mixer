# UGMG and BGMG - Univariate and Bivariate Gaussian Mixtures for GWAS summary statistics - `v1.0`

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
wget https://www.dropbox.com/s/5cvqzaayg1tn5u0/all_chromosomes_multiallelic_replaced_with_AT.ref.gz?dl=1 -O all_chromosomes_multiallelic_replaced_with_AT.ref.gz
wget https://www.dropbox.com/s/lxjwc5ub4yblqz0/BGMG_reference_data_11Msnps_2018_02_12.mat?dl=1 -O BGMG_reference_data_11Msnps_2018_02_12.mat
wget https://www.dropbox.com/s/j08f848raabcmcf/all_chromosomes.ref.gz?dl=1 -O all_chromosomes.ref.gz
```
``BGMG_reference_data_11Msnps.mat`` contains information about linkage disequilibrium structure, extracted for EUR population from 1000 genomes project, phase 3. Additional details are given further below.
``all_chromosomes_multiallelic_replaced_with_AT.ref.gz`` contains a list of 11015833 variants that we use to fit UGMG and BGMG model, with additional hack to replace multiallelic SNPs with ``A/T`` to workaround limitations of ``sumstats.py`` script. The original reference file is also enclosed (all_chromosomes.ref.gz). It was produced after basic QC steps, applied to all SNPs from 1000 genomes project, phase 3. 

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
matlab -nodisplay -nosplash -nodesktop -r "trait1='PGC_SCZ_2014.mat'; reference_data='BGMG_reference_data_11Msnps_2018_02_12.mat'; data_path='.'; out_file='PGC_SCZ_2014.result'; BGMG_run; exit;"
```
Alternatively, you may open matlab, create variables ``trait1``, ``reference_data``, ``data_path``, ``out_file`` with values as in the script above, and then execute ``BGMG_run.m`` script. The meaning of the parameters is as follows:
* ``trait1`` points to the output of ``sumstats.py``, e.i. summary stats converted to matlab format. It must be a filename, without path, of the ``.mat`` file containing ``zvec`` of the trait to be analyzed.
* ``data_path`` absolute or relative location of ``trait1`` file. Set this to ``.`` if file is located in current working directory.
* ``reference_data`` should be set to absolute or relative path of ``BGMG_reference_data_11Msnps_2018_02_12.mat``
* ``out_file`` specifies the name of the resulting file

## Interpret univariate results (UGMG)

After execution you should see approximately the following log file,
and a ``.pdf`` file with the following [QQ plot](https://github.com/precimed/BGMG/blob/master/PGC_SCZ_2014.png).

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
```

## Running bivariate model (BGMG)

TBD.

## Interpret bivariate results (BGMG)

TBD.

## ``BGMG_reference_data_11Msnps.mat`` file

TBD - add more details.

``BGMG`` comes with pre-calculated LD scores for european populations, stored in ``BGMG_reference_data_11Msnps.mat``.
The file contains the following information:
* chrnumvec - chromosome label, 1 to 22, for each variant from the reference
* posvec - base-pair position, build hg19, one value for each variant
* mafvec - minor allele frequency
* tldvec - total LD score
* numSNPsInLDr2_gt_r2min_vec - size of LD block
* hapmap3_mask - boolean mask of hapmap3 SNPs (to constrain inference to 1.1M SNPs, same as in LD score regression)
* LDr2_p8sparse - binary LD matrix for random pruning, thresholded at ``r2 >= 0.8``
* ref_ld - additional information about LD structure

The later field, ``ref_ld``, collected using ``https://github.com/ofrei/ldsc.git``, a modified version of LD Score Regression code,
to calculate LD scores (sum of squared allelic calculation, `r2`, and sum of its fourth power, `r4`).

## Authors

* Dominic Holland (Center for Multimodal Imaging and Genetics, University of California at San Diego)
* Oleksandr Frei (NORMENT, University of Oslo)
* Anders M. Dale (Center for Multimodal Imaging and Genetics, University of California at San Diego)
