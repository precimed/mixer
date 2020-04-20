MiXeR v1.2
==========

* MATLAB code is removed, from now you should run MiXeR via python interface. This is now described in [README](README.md) in the root of MiXeR github repository.
* MiXeR no longer provides standard errors on bivariate parameter estimates.
  Univariate parameter estimates can still, in principle, be calculated with ``--ci-alpha 0.05`` option, but this is turned off default, and is not recommended.
  Inference about statistical significance of the parameters should be made in context of the model selection criterion: Akaike (AIC) or Bayesian (BIC).
* You are required to update your scripts and use new reference files, as described below. 
  There is no need to update input files with summary statistics.
* Previous MiXeR version is available from git tag [v1.1](https://github.com/precimed/mixer/tree/v1.1).
  ``git checkout tags/v1.1 -b master`` allows to fetch it, more info [here](https://stackoverflow.com/questions/791959/download-a-specific-tag-with-git).
* After upgrading to the new version of the MiXeR, remember to re-compile native MiXeR plugin (``libbgmg.so``). Build instructions are available in the [README](README.md).
* ``mixer.py fit`` is split into four commands: ``fit1``, ``test1`` (univariate analysis) and ``fit2``, ``test2`` (bivariate analysis).
  ``fit`` commands are time-consuming and estimate MiXeR model parameters, typically from a subset of SNPs such as HapMap3 (can be passed with ``--extract`` flag).
  ``test`` commands are much faster and estimate power plots, QQ-plots, and similar, typically on a full set of SNPs.
  All commands have custom defaults, so you no longer have to specify ``--fit-sequence``, ``--kmax``, ``--qq-plots``, ``--power-curves``, etc.
  See the [README](README.md) file for examples.
* MiXeR estimates in bivariate analysis may, in some cases, depend on which of the two traits is selected as ``--trait1``.
  Therefore, it may be reasonable to run MiXeR twice, swapping the traits. However, the main reason for such instability 
  is due to low power in GWAS, and AIC/BIC will also capture such an issue. Therefore, if you observe significant
  AIC/BIC values, it may not be needed to run the MiXeR model twice.
* ``--plink-ld-bin0`` and ``--frq-file`` arguments are removed. Use ``--ld-file`` argument instead. 
  You need to download a new reference files for these arguments, links provided in the [README](README.md).
  Internally, the new references file are available on NIRD ``<SUMSTAT>/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld``.
  Further, there is an internal reference based on UK Biobank.
* ``mixer_figures.py`` now generates more data, including AIC and BIC, and also Dice coefficient for bivariate analysis.
  It also generates 
* Added support for Intel compiler, with build instructions for SAGA
* Added ``mixer.py perf`` option to evaluate performance of MiXeR cost funciton 
* ``mixer.py ld`` has simpler procedure to generate files for ``--ld-file`` argument based on your custom genotype panel 
  (plink is no longer required).
* ``--z1max`` and ``--z2max``allow to specify thresholds for right-censoring.
* MiXeR v1.2 uses precicely the same fit procedure as MiXeR v1.1, however it has a certain changes in how the LD matrix is stored internally, and also in how the cost function is evaluated. Particularly, the bivariate cost function based on sampling approach is now stateless, and it now consuming substentially less memory. The univariate cost function is based on convolution approach, as in v1.1, however it is now coupled with sampling for large z-scores to allow for ``--z1max`` and ``--z2max`` parameters.
  
