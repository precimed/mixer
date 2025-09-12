# MiXeR on simulated data

This usecase describe how to run MiXeR analysis (http://github.com/precimed/mixer) on synthetic data, generated with simu_linux (http://github.com/precimed/simu).

Prerequisites:

* have a local copy of https://github.com/comorment/mixer repository; below this folder is referred to as ``$MIXER_REFERENCE_FOLDER``
* have a local copy of mixer.sif container; below full path to this container is referred to as ``$MIXER_SIF``

Step 0. Setup ``${SINGULARITY_BIND}`` as follows:
```
export SINGULARITY_BIND=${MIXER_REFERENCE_FOLDER}/reference:/REF
```

Note that you may bypass the first two steps by using six ``.sumstats.gz`` files from ([this folder](https://github.com/comorment/mixer/tree/main/reference/hapgen)).
If you do this, go straight to step 3.

Step 1. Generate two pairs of traits for MiXeR analysis, covering two scenarios: complete polygenic overlap (``shared``), non-overlaping causal variants (``unique``), and partly overlaping causal variants (``partial``):
```
singularity shell --home $PWD:/home $SIF/mixer.sif 
simu_linux --qt --bfile /REF/hapgen/chr21 --seed 123 --causal-n 100 100 --trait1-sigsq 1 0 --trait2-sigsq 0 1 --num-components 2 --out unique --num-traits 2
simu_linux --qt --bfile /REF/hapgen/chr21 --seed 123 --causal-n 100     --trait1-sigsq 1   --trait2-sigsq 1   --num-components 1 --out shared --num-traits 2
simu_linux --qt --bfile /REF/hapgen/chr21 --seed 123 --causal-n 50 50 50 --trait1-sigsq 1 1 0   --trait2-sigsq 1 0 1   --num-components 3 --out partial --num-traits 2
plink2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno shared.pheno --pheno-name trait1 --out shared.1
plink2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno shared.pheno --pheno-name trait2 --out shared.2
plink2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno unique.pheno --pheno-name trait1 --out unique.1
plink2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno unique.pheno --pheno-name trait2 --out unique.2
plink2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno partial.pheno --pheno-name trait1 --out partial.1
plink2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno partial.pheno --pheno-name trait2 --out partial.2
```

Step 2. Convert the output of ``plink2`` to a format compatible with ``MiXeR``:
```
singularity exec --home $PWD:/home $SIF/mixer.sif python
import pandas as pd
for fname, out in [('{x}.{y}.trait{y}.glm.linear'.format(x=x,y=y), '{}.{}.sumstats'.format(x,y)) for x in ['unique', 'shared', 'partial'] for y in ['1', '2']]:
    pd.read_csv(fname, delim_whitespace=True)[['ID', 'REF', 'ALT', 'OBS_CT', 'T_STAT']].rename(columns={'ID':'SNP', 'REF':'A1', 'ALT':'A2', 'OBS_CT':'N', 'T_STAT':'Z'}).to_csv(out,index=False, sep='\t')
```

Step 3. Run univariate MiXeR analysis (each analysis should take up to 10 minutes on a standard laptop, in this demo example)

```
singularity shell --home $PWD:/home $SIF/mixer.sif

export MIXER_COMMON_ARGS="--chr2use 21 --ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep1.snps"

python /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS --trait1-file unique.1.sumstats --out unique.1

# now you have two options
# Option 1. Run the above command for the other 5 traits:
python /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS --trait1-file unique.2.sumstats --out unique.2
python /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS --trait1-file shared.1.sumstats --out shared.1
python /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS --trait1-file shared.2.sumstats --out shared.2
python /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS --trait1-file partial.1.sumstats --out partial.1
python /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS --trait1-file partial.2.sumstats --out partial.2

# Option 2. Do a bit of cheating, and copy copy the .json files, enforcing exacly the same genetic architecture 
# ("pi", "sig2beta" and "sig2zero" parameters) in all six traits. This shouldn't make much difference as long as 
# univariate estimates are fairly accurate.
cp unique.1.json unique.2.json && cp unique.1.json shared.1.json && cp unique.1.json shared.2.json && cp unique.1.json partial.1.json && cp unique.1.json partial.2.json

# Now proceed with bivariate analyses
python /tools/mixer/precimed/mixer.py fit2 $MIXER_COMMON_ARGS --trait1-file unique.1.sumstats --trait2-file unique.2.sumstats --trait1-params unique.1.json --trait2-params unique.2.json --out unique

python /tools/mixer/precimed/mixer.py fit2 $MIXER_COMMON_ARGS --trait1-file shared.1.sumstats --trait2-file shared.2.sumstats --trait1-params shared.1.json --trait2-params shared.2.json --out shared

python /tools/mixer/precimed/mixer.py fit2 $MIXER_COMMON_ARGS --trait1-file partial.1.sumstats --trait2-file partial.2.sumstats --trait1-params partial.1.json --trait2-params partial.2.json --out partial

python /tools/mixer/precimed/mixer_figures.py two --json shared.json --out shared
python /tools/mixer/precimed/mixer_figures.py two --json unique.json --out unique
python /tools/mixer/precimed/mixer_figures.py two --json partial.json --out partial
```

Step 4. The results are available in ``shared.png``, ``unique.png`` and ``partial.png`` figures.

Step 5. You could also use [MIXER_SIMU.job](mixer_simu/MIXER_SIMU.job) script to run this on a cluster 20 times, 
each with a random subsets of SNPs as defined by ``1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep[1-20].snps`` files, 
and use variation in parameter estimates amoung these runs this to infer uncertainty of parameter estimates.
Adjust ``MIXER_SIMU.job`` to match configuration of your cluster. Then trigger analysis by running ``sbatch MIXER_SIMU.job`` command.
Once you get all results, combine them across 20 runs, and produce figures as follows:

```
singularity shell --home $PWD:/home $SIF/mixer.sif
python /tools/mixer/precimed/mixer_figures.py combine --json unique.fit.rep@.json --out unique.fit
python /tools/mixer/precimed/mixer_figures.py combine --json unique.test.rep@.json --out unique.test
python /tools/mixer/precimed/mixer_figures.py combine --json shared.fit.rep@.json --out shared.fit
python /tools/mixer/precimed/mixer_figures.py combine --json shared.test.rep@.json --out shared.test
python /tools/mixer/precimed/mixer_figures.py combine --json partial.fit.rep@.json --out partial.fit
python /tools/mixer/precimed/mixer_figures.py combine --json partial.test.rep@.json --out partial.test

python /tools/mixer/precimed/mixer_figures.py two --statistic mean std --json-fit unique.fit.json --json-test unique.test.json --out mixer_simu_unique
python /tools/mixer/precimed/mixer_figures.py two --statistic mean std --json-fit shared.fit.json --json-test shared.test.json --out mixer_simu_shared
python /tools/mixer/precimed/mixer_figures.py two --statistic mean std --json-fit partial.fit.json --json-test partial.test.json --out mixer_simu_partial
```
After processing the resulting figures shoud look like this:

Unique:
![mixer_simu_unique.png](https://raw.githubusercontent.com/precimed/mixer/master/usecases/mixer_simu/mixer_simu_unique.png)

Partial:
![mixer_simu_partial.png](https://raw.githubusercontent.com/precimed/mixer/master/usecases/mixer_simu/mixer_simu_partial.png)

Shared:
![mixer_simu_shared.png](https://raw.githubusercontent.com/precimed/mixer/master/usecases/mixer_simu/mixer_simu_shared.png)

To include density plots to the resulting figure, add ``--trait1-file`` and ``--trait2-file`` arguments like this:
```
python /tools/mixer/precimed/mixer_figures.py  two --statistic mean std --json-fit unique.fit.json --json-test unique.test.json  --out mixer_simu_unique --trait1-file unique.1.sumstats --trait2-file unique.2.sumstats
```

## Docker details

The syntax above can be adapted to Docker as follows:

```bash
#!/bin/bash
# define environment variables:
export MIXER_DOCKER_IMAGE="ghcr.io/precimed/gsa-mixer:latest"  # adapt as necessary
export MIXER_REFERENCE_FOLDER=<full path> # adapt as necessary
export SLURM_ARRAY_TASK_ID=1
export MIXER_COMMON_ARGS="--chr2use 21 --ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim"
export REP="rep${SLURM_ARRAY_TASK_ID}"
export EXTRACT="--extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.$REP.snps"
# shortcuts for different binaries and interactive shell:
export PYTHON="docker run --platform=linux/amd64 --rm -v ${PWD}:/home -v ${MIXER_REFERENCE_FOLDER}/reference:/REF -w/home --entrypoint=python ${MIXER_DOCKER_IMAGE}"
export SIMU="docker run --platform=linux/amd64 --rm -v ${PWD}:/home -v ${MIXER_REFERENCE_FOLDER}/reference:/REF -w/home --entrypoint=simu_linux ${MIXER_DOCKER_IMAGE}"
export PLINK1="docker run --platform=linux/amd64 --rm -v ${PWD}:/home -v ${MIXER_REFERENCE_FOLDER}/reference:/REF -w/home --entrypoint=plink1 ${MIXER_DOCKER_IMAGE}"
export PLINK2="docker run --platform=linux/amd64 --rm -v ${PWD}:/home -v ${MIXER_REFERENCE_FOLDER}/reference:/REF -w/home --entrypoint=plink2 ${MIXER_DOCKER_IMAGE}"
export ISHELL="docker run --platform=linux/amd64 --rm -it -v ${PWD}:/home -v ${MIXER_REFERENCE_FOLDER}/reference:/REF -w/home --entrypoint=bash ${MIXER_DOCKER_IMAGE}"

# invoke mixer.py help documentation:
$PYTHON /tools/mixer/precimed/mixer.py fit1 --help 

# prepare simulated data:
$SIMU --qt --bfile /REF/hapgen/chr21 --seed 123 --causal-n 100 100 --trait1-sigsq 1 0 --trait2-sigsq 0 1 --num-components 2 --out unique --num-traits 2
$SIMU --qt --bfile /REF/hapgen/chr21 --seed 123 --causal-n 100     --trait1-sigsq 1   --trait2-sigsq 1   --num-components 1 --out shared --num-traits 2
$SIMU --qt --bfile /REF/hapgen/chr21 --seed 123 --causal-n 50 50 50 --trait1-sigsq 1 1 0   --trait2-sigsq 1 0 1   --num-components 3 --out partial --num-traits 2
$PLINK2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno shared.pheno --pheno-name trait1 --out shared.1
$PLINK2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno shared.pheno --pheno-name trait2 --out shared.2
$PLINK2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno unique.pheno --pheno-name trait1 --out unique.1
$PLINK2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno unique.pheno --pheno-name trait2 --out unique.2
$PLINK2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno partial.pheno --pheno-name trait1 --out partial.1
$PLINK2 --bfile /REF/hapgen/chr21 --glm allow-no-covars --pheno partial.pheno --pheno-name trait2 --out partial.2

$PYTHON -c "import pandas as pd
for fname, out in [('{x}.{y}.trait{y}.glm.linear'.format(x=x,y=y), '{}.{}.sumstats'.format(x,y)) for x in ['unique', 'shared', 'partial'] for y in ['1', '2']]:
    pd.read_csv(fname, delim_whitespace=True)[['ID', 'REF', 'ALT', 'OBS_CT', 'T_STAT']].rename(columns={'ID':'SNP', 'REF':'A1', 'ALT':'A2', 'OBS_CT':'N', 'T_STAT':'Z'}).to_csv(out,index=False, sep='\t')"

# MiXeR univariate fit:
$PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file unique.1.sumstats --out unique.1

# cheat:
cp unique.1.json unique.2.json && cp unique.1.json shared.1.json && cp unique.1.json shared.2.json && cp unique.1.json partial.1.json && cp unique.1.json partial.2.json

# bivariate fit:
$PYTHON /tools/mixer/precimed/mixer.py fit2 $MIXER_COMMON_ARGS $EXTRACT --trait1-file unique.1.sumstats --trait2-file unique.2.sumstats --trait1-params unique.1.json --trait2-params unique.2.json --out unique
$PYTHON /tools/mixer/precimed/mixer.py fit2 $MIXER_COMMON_ARGS $EXTRACT --trait1-file shared.1.sumstats --trait2-file shared.2.sumstats --trait1-params shared.1.json --trait2-params shared.2.json --out shared
$PYTHON /tools/mixer/precimed/mixer.py fit2 $MIXER_COMMON_ARGS $EXTRACT --trait1-file partial.1.sumstats --trait2-file partial.2.sumstats --trait1-params partial.1.json --trait2-params partial.2.json --out partial

# figures:
$PYTHON /tools/mixer/precimed/mixer_figures.py two --json shared.json --out shared
$PYTHON /tools/mixer/precimed/mixer_figures.py two --json unique.json --out unique
$PYTHON /tools/mixer/precimed/mixer_figures.py two --json partial.json --out partial
```
