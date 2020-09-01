MiXeR on TSD
============

- load "module load CMake/3.15.3-GCCcore-8.3.0 Boost/1.73.0-GCCcore-8.3.0 Python/3.7.4-GCCcore-8.3.0"
- download https://github.com/precimed/mixer/archive/master.zip, import to TSD, place on /cluster and extract
- mkdir <MIXER_ROOT>/src/build && cd <MIXER_TOOL>/src/build && cmake .. && make -j8
- adjust <MIXER_ROOT>/scripts/tsd_ugmg_script.sh and <MIXER_ROOT>/scripts/tsd_bgmg_script.sh   - not fully tested, minor changes may be needed here. Submit with "sbatch <script.sh>" as usual, it will internally trigger a job array with 20 runs.
