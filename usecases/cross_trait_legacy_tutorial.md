# Legacy tutorial on cross-trait MiXeR

## Contents

* [Introduction](#introduction)
* [Data downloads](#data-downloads)
* [Data preparation](#data-preparation)
* [Run MiXeR](#run-mixer)
* [MiXeR options](#mixer-options)
* [Visualize MiXeR results](#visualize-mixer-results)

## Introduction

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
  * We recommed to use our own scripts to pre-process summary statistcs (clone from [here](https://github.com/precimed/python_convert)):
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
  * [DEPRECATED] If you use MiXeR v1.1 and v1.2, it will recognize summary statistics in LDSC format as described [here](https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format). In brief, each trait must be represented as a single table containing columns SNP, N, Z, A1, A2. Thus, it is possible to use ``munge_sumstats.py`` script as described [here](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability#step-1-download-the-data). This might be convenient for users who are already familiar with LDSR functionality. In MiXeR v1.3 the format for GWAS summary statistics files is the same, but now I advice against using HapMap3 to constrain the set of SNPs for the fit procedure. Instead, MiXeR should receive a full set of summary statistics from a GWAS on imputed genotype data. Also, note that for case/control ``munge_sumstats.py`` from LD Score Regression generate sample size as a sum ``n = ncase + ncontrol``. We recommend to use ``neff = 4 / (1/ncase + 1/ncontrol)`` to account for imbalanced classes.
  
* Generate ``.ld`` and ``.snps`` files from the reference panel. Note that this step optional if you work with EUR-based summary statistics. To use EUR reference, simply download the files from  [here](https://github.com/comorment/mixer/tree/main/reference/ldsc/1000G_EUR_Phase3_plink) or take it from NIRD (``/projects/NS9114K/MMIL/SUMSTAT/LDSR/1000G_EUR_Phase3_plink``) or TSD (`` ``) if you have access. NB! Download size is around ``24 GB``.
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
  
  * To generate ``1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.repNN.snps`` files, repeat the following in a loop for
    ```
    export REP=1 # repeat for REP in 1..20
    python3 <MIXER_ROOT>/precimed/mixer.py snps \
       --lib <MIXER_ROOT>/src/build/lib/libbgmg.so \
       --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
       --ld-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
       --out LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps \
       --maf 0.05 --subset 2000000 --r2 0.8 --seed ${SLURM_ARRAY_TASK_ID}
    ```
    This can be done as a job array.
    Note that in the code above the ``@`` symbol does NOT need to be replace with an actual chromosome. It should stay as ``@`` in your command.
    
## Run MiXeR

### Univariate analysis

Fit the model:
```
python3 <MIXER_ROOT>/precimed/mixer.py fit1 \
      --trait1-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --out SSGAC_EDU_2018_no23andMe_noMHC.fit.rep${SLURM_ARRAY_TASK_ID} \
      --extract LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  <MIXER_ROOT>/src/build/lib/libbgmg.so \
```

Apply the model to the entire set of SNPs, without constraining to ``LDSR/w_hm3.justrs``:
```
python3 <MIXER_ROOT>/precimed/mixer.py test1 \
      --trait1-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --load-params-file SSGAC_EDU_2018_no23andMe_noMHC.fit.rep${SLURM_ARRAY_TASK_ID}.json \
      --out SSGAC_EDU_2018_no23andMe_noMHC.test.rep${SLURM_ARRAY_TASK_ID} \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  <MIXER_ROOT>/src/build/lib/libbgmg.so \
```

The results will be saved ``<out_file>.json`` file.
Repeat the above analysis for the second trait (``PGC_SCZ_2014_EUR_qc_noMHC.csv.gz``).

To visualize the results:
```
python precimed/mixer_figures.py combine --json PGC_SCZ_2014_EUR_qc_noMHC.fit.rep@.json --out combined/PGC_SCZ_2014_EUR.fit
python precimed/mixer_figures.py one --json PGC_SCZ_2014_EUR.json --out PGC_SCZ_2014_EUR --statistic mean std
```

### Bivariate (cross-trait) analysis

Fit the model:
```
python3 <MIXER_ROOT>/python/mixer.py fit2 \
      --trait1-file PGC_SCZ_2014_EUR_qc_noMHC.csv.gz \
      --trait2-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --trait1-params-file PGC_SCZ_2014_EUR_qc_noMHC.fit.rep${SLURM_ARRAY_TASK_ID}.json \
      --trait2-params-file SSGAC_EDU_2018_no23andMe_noMHC.fit.rep${SLURM_ARRAY_TASK_ID}.json \
      --out PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.fit.rep${SLURM_ARRAY_TASK_ID} \
      --extract LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  <MIXER_ROOT>/src/build/lib/libbgmg.so \
```

Apply the model to the entire set of SNPs, without constraining to ``LDSR/w_hm3.justrs``:
```
python3 <MIXER_ROOT>/python/mixer.py test2 \
      --trait1-file PGC_SCZ_2014_EUR_qc_noMHC.csv.gz \
      --trait2-file SSGAC_EDU_2018_no23andMe_noMHC.csv.gz \
      --load-params-file PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.fit.rep${SLURM_ARRAY_TASK_ID}.json \
      --out PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC.test.rep${SLURM_ARRAY_TASK_ID} \
      --bim-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  <MIXER_ROOT>/src/build/lib/libbgmg.so \
```

Note that these parameters point to the results of univariate analysis for both traits, so those must be generated first.
The results will be saved ``<out_file>.json`` file.

To visualize the results (where `<prefix>` is, for example, ``PGC_SCZ_2014_EUR_qc_noMHC_vs_SSGAC_EDU_2018_no23andMe_noMHC``):
```
python precimed/mixer_figures.py combine --json <prefix>.fit.rep@.json --out <prefix>.fit
python precimed/mixer_figures.py combine --json <prefix>.test.rep@.json --out <prefix>.test
python precimed/mixer_figures.py two --json-fit <prefix>.fit.json --json-test <prefix>.test.json --out <out_file> --statistic mean std
```

## MiXeR options

Run ``--help`` commands to list available options and their description.
```
python3 mixer.py ld --help
python3 mixer.py snps --help
python3 mixer.py fit1 --help
python3 mixer.py test1 --help
python3 mixer.py fit2 --help
python3 mixer.py test2 --help
python3 mixer.py perf --help
```

## Visualize MiXeR results

First step is to average the results across 20 runs with ``mixer_figures.py combine``, which works for univariate (fit1, test1) and bivariate (fit2, test2) runs:
```
python precimed/mixer_figures.py combine --json <prefix>.fit.rep@.json --out <prefix>.fit
```

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

These output is described in the cross-trait MiXeR publication. 
