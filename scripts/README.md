MiXeR on TSD
============

- ``module load CMake/3.15.3-GCCcore-8.3.0 Boost/1.73.0-GCCcore-8.3.0 Python/3.7.4-GCCcore-8.3.0``
- ``cd ~ && git clone --recurse-submodules -j8 https://github.com/precimed/mixer.git && tar -czvf mixer_master.tar.gz mixer ``, then import to TSD, place on /cluster and extract. Note that "Download" button on github won't work correctly, because it does not include submodules (i.e. zlib) in the package.
- ``mkdir <MIXER_ROOT>/src/build && cd <MIXER_TOOL>/src/build && cmake .. && make -j8``
- adjust ``<MIXER_ROOT>/scripts/tsd_ugmg_script.sh`` and ``<MIXER_ROOT>/scripts/tsd_bgmg_script.sh`` (not fully tested, minor changes may be needed here).
- Submit with "sbatch <script.sh>" as usual, it will internally trigger a job array with 20 runs.
- When results are ready, create figures as described in the main README.md file (see "visualize the results" section)
- first-time configuration of your Python environment:
  ```
  ssh p33-appn-norment01
  module load  Python/3.7.4-GCCcore-8.3.0

  # clean install jupyter in a new environment
  python3 -m venv /cluster/projects/p33/users/<user>/py3
  source /cluster/projects/p33/users/<user>/py3/bin/activate   # best to have this on /cluster

  pip3 install --index-url=file:///shared/pypi/mirror/web/simple numpy
  pip3 install --index-url=file:///shared/pypi/mirror/web/simple pandas
  pip3 install --index-url=file:///shared/pypi/mirror/web/simple scipy==1.2.0rc1
  pip3 install --index-url=file:///shared/pypi/mirror/web/simple matplotlib
  pip3 install --index-url=file:///shared/pypi/mirror/web/simple statsmodels
  pip3 install --index-url=file:///shared/pypi/mirror/web/simple numdifftools
  pip3 install --index-url=file:///shared/pypi/mirror/web/simple matplotlib_venn
  pip3 install --index-url=file:///shared/pypi/mirror/web/simple jupyter
  ```
  Then use this Python environment as follows:
  ```
  module load  Python/3.7.4-GCCcore-8.3.0
  source /cluster/projects/p33/users/<user>/py3/bin/activate
  ```

