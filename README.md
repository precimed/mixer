# BGMG (Bivariate Gaussian Mixture for GWAS summary statistics) `v1.0`

`BGMG` is a set of matlab scripts for estimating polygenicity and polygenic overlap between two traits given their GWAS summary statistics.

## Getting Started

To download `BGMG` you should clone these two repositories:
```  
git clone https://github.com/precimed/bgmg.git
git clone https://github.com/precimed/python_convert.git
```

``BGMG`` repository contains the actual matlab code of the model.
``python_convert`` repository contains supplementary code that loads summary statistics into ``BGMG``-compatible format.

In addition, you may want to clone ``https://github.com/ofrei/ldsc.git``, a modified version of LD Score Regression code,
to calculate LD scores (sum of squared allelic calculation, `r2`, and sum of its fourth power, `r4`).
If you summary statistics came from european population you may proceed using LD scores from ``reference_data`` folder, included in this repository. In this case  you DO NOT need ``ofrei/ldsc`` code, but you still need ``precimed/python_convert`` code to load summary statistics into BGMG format.

## Prepare summary statistics

First step is to load summary statistics into `BGMG`-compatible format, using ``sumstats.py`` tool from ``python_convert`` repository.
The following example shows how to load Schizophrenia and Triglycerides level GWAS from PGC and LIPIDS consortia:

```
python sumstats.py csv \
	--sumstats scz2.snp.results.txt.gz
	--ncase-val 35476 \
	--ncontrol-val 46839 \
	--auto  \
	--chr hg19chrc \
	--out PGC_SCZ_2014.csv && gzip PGC_SCZ_2014.csv \

python sumstats.py csv \
	--sumstats jointGwasMc_TG.txt.gz
	--auto  \
	--chrpos SNP_hg19 \
	--out LIPIDS_TG_2013.csv && gzip LIPIDS_TG_2013.csv \
  
python sumstats.py mat --sumstats PGC_SCZ_2014.csv.gz   --ref reference_data/1m.ref.gz --out PGC_SCZ_2014.mat
python sumstats.py mat --sumstats LIPIDS_TG_2013.csv.gz --ref reference_data/1m.ref.gz --out LIPIDS_TG_2013.mat
```
For more information about ``sumstats.py`` utility see ``python_convert`` repository,
or type ``sumstats.py --help`` or ``sumstats.py <command> --help`` where ``<command>`` is either ``csv`` or ``mat``.

This results in a separate ``.mat`` file for each trait.
Each file contains two variables, ``zvec`` and ``nvec`` ---
a signed test statistic (`z-score`) and effective sample size (`nvec`) for each variant. 
``zvec`` and ``nvec`` will be a column vector of length `N` (number of variants from the reference file, ``1m.ref.gz``).
The order of the variants will match the reference file.
``zvec`` or ``nvec`` might contain undefined values, ``nan``.
In this case the analysis will be limited to a subset of variants where ``zvec`` and ``nvec`` are defined across both traits.

## Running BGMG model

You may run ``BGMG`` from console as follows:
```
matlab -nodisplay -nosplash -nodesktop -r "trait1='PGC_SCZ_2014.mat'; trait2='LIPIDS_TG_2013.mat'; data_path='.'; BGMG_run; exit;"
```
Here ``data_path`` refers to location of ``PGC_SCZ_2014.mat`` and ``LIPIDS_TG_2013.mat`` files.
Set ``data_path=''`` if you choose to specify full path for ``trait1`` and ``trait2`` files.

Remember to set your working directory to the folder containing ``BGMG_run``.

Other supported options:

1. Use ``reference_data`` to specify path of your custom reference (see below for instructions how to generate one):
```
matlab -nodisplay -nosplash -nodesktop -r "trait1='PGC_SCZ_2014.mat'; trait2='LIPIDS_TG_2013.mat'; data_path='.'; reference_data='<full_path_to_reference_data_folder>'; BGMG_run; exit;"
```

## Interpret BGMG results

TBD.

Results of ``BGMG_run`` script are stored in a text file, named after input files, for example ``BGMG_run_PGC_SCZ_2014-LIPIDS_TG_2013.txt``. It includes univariate analysis for each of the traits, and bivariate cross-traits analysis.

## Estimating polygenicity via univariate analysis (UGMG model)

``BGMG_run`` script also supports univariate model, triggered as follows:

```
matlab -nodisplay -nosplash -nodesktop -r "trait1='PGC_SCZ_2014.mat';   BGMG_run; exit;"
matlab -nodisplay -nosplash -nodesktop -r "trait1='LIPIDS_TG_2013.mat'; BGMG_run; exit;"
```

The result is stored in files ``UGMG_run_PGC_SCZ_2014.txt`` and ``UGMG_run_LIPIDS_TG_2013.txt``.

## ``reference_data`` folder for EUR populations.

``BGMG`` comes with pre-calculated LD scores for european populations.
LD scores are calculated for ``1190321`` markers from ``reference_data\1m.ref.gz`` file,
defined as overlap between hapmap3 and 1kG phase3 markers.
LD scores are calculated towards all markers, available in 1kG phase3 data,
except for ``w_ld`` scores which are calculated withing ``1m.ref`` template (e.g. towareds ``1190321`` markers).

Input data used in LD score estimation:
* https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 --- list of hapmap3 SNPs
* https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz --- 1kG phase 3 genotypes, 489 EUR individuals (14 excluded due to relatedness)

The procedure is self-described by ``reference_data/Makefile``.
To run ``Makefile`` you need ``GNU Make`` utility.
Note that in this case ``make`` is not used to build ``C++`` code,
but rather as a fancy alternative to writting shell scripts --- just to automate the process.
You can also run it in parallel with ``make -j8``.

## Requirements

1. ``Matlab`` (tested with ``8.6.0.267246 (R2015b)``, may work with other versions as well)
2. ``python 3 > version >= 2.7`` --- to estimate 

## Authors

Oleksandr Frei (NORMENT, University of Oslo)

Anders M. Dale (Center for Multimodal Imaging and Genetics, University of California at San Diego)
