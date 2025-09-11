# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### [Added]

- Examples of how to apply MiXeR on [synthetic](usecases/mixer_real.md) and [real](usecases/mixer_simu.md) data
- Added [`usecases/run_mixer.ipynb`](usecases/run_mixer.ipynb) Jupyter notebook for simplifying Slurm job submissions for lists of traits.
- Include ``module purge`` in Slurm job scripts to avoid conflicts with system modules.
- Added example outputs from GSA-MiXeR analysis (`usecases/gsa_mixer/out_example`)
- Example scripts for GSA-MiXeR analysis (``scripts/``)

### [Fixed]

- Removed use of `np.unique` in [`usecases/run_mixer.ipynb`](usecases/run_mixer.ipynb)
- Fixed incorrect job dependency bug in [`usecases/run_mixer.ipynb`](usecases/run_mixer.ipynb)
- Updated outdated paths to comorment containers on p697
- Explicitly named log files in Slurm .job scripts
- Fixed reference to wrong container file in MIXER_SIMU.job

