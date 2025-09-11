# MiXeR

This project provides [Singularity](https://sylabs.io/singularity/) wrapper for <https://github.com/precimed/mixer>, 
as well as [Docker](https://www.docker.com/) images for the same software.

The original reference and example data are included in this repository.

If you use this package, please cite the original work and all reference data used in your analysis.

The history of changes is available in the [CHANGELOG.md](CHANGELOG.md) file.

To get started, see:

- general instructions about how to use CoMorMent containers: <https://github.com/comorment/containers#getting-started>

- examples of how to use mixer.sif container on [synthetic](usecases/mixer_simu.md) and [real](usecases/mixer_real.md) data.

- Jupyter notebook [usecases/run_mixer.ipynb](usecases/run_mixer.ipynb) for submitting a bunch of MiXeR jobs at once, with lists of primary and secondary traits. 

## Installation and set up

### Dependencies on host system

In order to set up these resource, some software may be required

- [Singularity/SingularityCE](https://sylabs.io/singularity/) or [Apptainer](https://apptainer.org)
- [Git](https://git-scm.com/)
- [Git LFS](https://git-lfs.com)
- [ORAS CLI](https://oras.land)

### Clone the repository

To download the last revision of this project, issue:

```bash
cd path/to/repositories
git clone --depth 1 https://github.com/comorment/mixer.git
cd mixer
git lfs pull  # pull "large" files
```

### Update the `mixer.sif` container

To obtain updated versions of the Singularity Image Format (.sif) container file `, issue

```bash
cd path/to/repositories/mixer/singularity
mv mixer.sif mixer.sif.old  # optional, just rename the old(er) file
apptainer pull docker://ghcr.io/comorment/mixer:<tag>  # or
singularity pull docker://ghcr.io/comorment/mixer:<tag> # or 
oras pull ghcr.io/comorment/mixer_sif:<tag>
```

where `<tag>` corresponds to a tag listed under [packages](https://github.com/comorment/mixer/pkgs/container/mixer), 
such as `latest`, `main`, or `sha_<GIT SHA>`. 
The `oras pull` statement pulls the `mixer.sif` file from [ghcr.io](https://github.com/comorment/mixer/pkgs/container/mixer_sif) using the [ORAS](https://oras.land) registry, without the need to build the container locally.

### Pulling and using Docker image

To pull the corresponding Docker image, issue:

```bash
docker pull ghcr.io/comorment/mixer:<tag>
```

If working on recent Macs, add the `--platform=linux/amd64` after `docker pull`. 
This may allow replacing `singularity exec ...` or `apptainer exec ...` statements with appropriate `docker run ...` statements in the [usecases/mixer_simu](usecases/mixer_simu#docker-details) section, 
on systems where Singularity or Apptainer is unavailable. 
Functionally, the Docker image is equivalent to the Singularity container, but note that syntax for mounting volumes and invoking commands may differ.
Please refer to[docs.docker.com](https://docs.docker.com) for more information.

> [!NOTE] Note that the provided Docker image may not support all CPUs, and may not be able to run on all systems via CPU virtualization.
> An option may be to build the Docker image on the host machine (e.g., M1 Macs, older Intel CPUs), as:
>
>```bash
>docker build --platform=linux/amd64 -t ghcr.io/comorment/mixer -f dockerfiles/mixer/Dockerfile .
>``` 

An exampe of using the Docker image is provided in the [usecases/mixer_simu](usecases/mixer_simu#docker-details) section.

## Systems without internet access

Some secure platforms do not have direct internet access, hence we recommend cloning/pulling all required files on a machine with internet access as explained above, and archive the `mixer` directory with all files and moving it using whatever file uploader is available for the platform.

```bash
cd /path/to/mixer
SHA=$(git rev-parse --short HEAD)
cd ..
tar --exclude=".git/*" -cvf mixer_$SHA.tar mixer
```

## Citation info

If you use the software provided here, please cite our relevant preprint:

```
Akdeniz, B.C., Frei, O., Hagen, E., Filiz, T.T., Karthikeyan, S., Pasman, J.A., Jangmo, A., Bergsted, J., Shorter, J.R., Zetterberg, R., Meijsen, J.J., Sønderby, I.E., Buil, A., Tesli, M., Lu, Y., Sullivan, P., Andreassen, O.A., & Hovig, E. (2022). COGEDAP: A COmprehensive GEnomic Data Analysis Platform. arXiv:2212.14103 [q-bio.GN]. DOI: [10.48550/arXiv.2212.14103](https://doi.org/)
```

Bibtex format:
```
@misc{akdeniz2022cogedap,
      title={COGEDAP: A COmprehensive GEnomic Data Analysis Platform}, 
      author={Bayram Cevdet Akdeniz and Oleksandr Frei and Espen Hagen and Tahir Tekin Filiz and Sandeep Karthikeyan and Joelle Pasman and Andreas Jangmo and Jacob Bergsted and John R. Shorter and Richard Zetterberg and Joeri Meijsen and Ida Elken Sonderby and Alfonso Buil and Martin Tesli and Yi Lu and Patrick Sullivan and Ole Andreassen and Eivind Hovig},
      year={2022},
      eprint={2212.14103},
      archivePrefix={arXiv},
      primaryClass={q-bio.GN}
}
```

Note that this project will soon fall under the "[COSGAP](https://cosgap.readthedocs.io/en/latest/)" umbrella, and that the citation info will be updated accordingly.

For the [MiXeR software](https://github.com/precimed/mixer) itself, if you use MiXeR software for your research publication, please cite the following paper(s):

* for univariate analysis: D. Holland et al., Beyond SNP Heritability: Polygenicity and Discoverability Estimated for Multiple Phenotypes with a Univariate Gaussian Mixture Model, PLOS Genetics, 2020, https://doi.org/10.1371/journal.pgen.1008612
* for cross-trait analysis: O.Frei et al., Bivariate causal mixture model quantifies polygenic overlap between complex traits beyond genetic correlation, Nature Communications, 2019, https://www.nature.com/articles/s41467-019-10310-0
