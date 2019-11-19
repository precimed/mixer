Cost function for Univariate and Bivariate mixture model
Written as a plain C interface. Implemented in C++. Require Boost libraries.

**Clone or Pull this repo**

This repository uses git submodules. To pull this repository including submodules se 
```
git clone --recurse-submodules -j8 git://github.com/precimed/bgmg.git
```
For more information about git submodules see
https://stackoverflow.com/questions/3796927/how-to-git-clone-including-submodules
https://stackoverflow.com/questions/1030169/easy-way-to-pull-latest-of-all-git-submodules

**Build on Windows**
- install MS VS 2015
- install pre-built boost libraries from https://sourceforge.net/projects/boost/files/boost-binaries/
- set env variables: 
```
setx BOOST_ROOT H:\local\boost_1_62_0
setx BOOST_LIBRARYDIR H:\local\boost_1_62_0\lib64-msvc-14.0
```
- ``cmake .. -G"Visual Studio 14 2015 Win64"``
- open bgmg.sln in Visual Studio, and build as usual

**Build on Linux**

It is important to build with reasonably  old version of gcc and boost libraries, otherwise it may conflict with MATLAB's (custom) version of C++ runtime (libstdc++.so.6 and libgcc_s.so.1). On Abel (UiO supercomputer) the following combination of modules turned out to work well:

```
module purge
module load plink
module load python2/2.7.10.gnu
module load matlab/R2017a
module load cmake
module load boost/1.55.0
```

Build:
``
mkdir build && cd build && cmake .. && make -j12
``

To run from matlab:
```
 trait1_file='/work/users/oleksanf/LDSR/LDSR_Data/PGC_SCZ_2014_noMHC.sumstats.gz';
 bim_file='/work/users/oleksanf/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim';
 frq_file='/work/users/oleksanf/LDSR/1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.@.frq';
 plink_ld_bin='/work/users/oleksanf/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin'; chr_labels = 1:22;
 out_file='/usit/abel/u1/dennisva/BGMG_results/PGC_SCZ_2014';
 bgmg_shared_library='src/build/lib/libbgmg.so'; bgmg_shared_library_header='src/bgmg_matlab.h';
 kmax=5000; max_causal_fraction=0.03; cache_tag_r2sum=0; SEED=123; randprune_r2=0.1; randprune_n=64; CI_ALPHA=0.05; r2min=0.0;
 DO_FIT_UGMG=1; FIT_FULL_MODEL=1; POWER_PLOT=1; POWER_PLOT_DOWNSCALE=100; QQ_PLOT=1; QQ_PLOT_DOWNSCALE=100; QQ_PLOT_BINS=1; QQ_PLOT_BINS_DOWNSCALE=50;
 UGMG_cpp_run_simple; 
```


Most important is to use boost/1.55.0, which cames with gcc/4.9.0. This will be compatible with matlab/R2017a which uses boost/1.56. 
Useful tricks to investigate versions:

```
ldd /cluster/software/VERSIONS/matlab/R2017a/bin/glnxa64/MATLAB
export "LD_LIBRARY_PATH=/cluster/software/VERSIONS/matlab/R2017a/sys/os/glnxa64:$LD_LIBRARY_PATH"
ldd libbgmg.so
```

Errors like ``cluster/software/VERSIONS/matlab/R2017a/sys/os/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found`` indicate that you've built libbgmg.so with too new version of gcc.

At MMIL servers we have quite old version of matlab ``matlab R2015b (8.6.0.267246)``, which cames with ``boost_1_49_0``, and ``GLIBCXX`` up to ``3.4.17``. This imply that we have to compile bgmglib.so with ``gcc 4.7.2`` or older. This compiler is not available at MMIL, so best way is to build this binary on Abel, and rsync to MMIL. So, the process is as follows:

1. ``module purge && module load cmake && module load gcc/4.7.2``
2. Download [boost 1.49.0](https://www.boost.org/doc/libs/1_49_0/more/getting_started/unix-variants.html) to abel (because this version is not available)
3. Build boost: ``./bootstrap.sh --with-libraries=program_options,filesystem,system,date_time && ./b2 --j12``
4. ``cd ~/BGMG_mmil/src/build && cmake .. -DBOOST_ROOT=/usit/abel/u1/oleksanf/boost_1_49_0``
5. At mmil: ``rsync -avzP oleksanf@abel.uio.no:/usit/abel/u1/oleksanf/precimed/BGMG_mmil/src/build/lib /home/oleksandr/precimed/BGMG/src/build``

UPDATE: for some reason I have to also rsync boost libraries, and put  them in LD_LIBRARY_PATH before starting matlba.
```
rsync -avzP oleksanf@abel.uio.no:/usit/abel/u1/oleksanf/boost_1_49_0/stage/lib
export "LD_LIBRARY_PATH=/home/oleksandr/precimed/BGMG/src/build/lib:$LD_LIBRARY_PATH"
```


**Test that build succeeded**

The build produces two important artifact: ``<build>/bin/bgmg-cli`` (executable) and ``<build>/lib/libbgmg.so`` (shared library).

To test ``bgmg-cli``:

```
./bgmg-cli --help
./bgmg-cli --bim /work/users/oleksanf/bfile_merged/chr@.bim --frq /work/users/oleksanf/bfile_merged/chr@.frq --trait1 XXX_CRP_2009_noMHC.sumstats.gz  & tailf bgmg.bgmglib.log
./bgmg-cli --bim /work/users/oleksanf/bfile_merged/chr@.bim --plink-ld bfile_merged_ldmat_p01_SNPwind50k_chr21.ld.gz --out bfile_merged_ldmat_p01_SNPwind50k_chr21.ld.bin & tailf bfile_merged_ldmat_p01_SNPwind50k_chr21.ld.bin.bgmglib.log
```

To test ``bgmglib.so`` it may be loaded to matlab:

```
bgmg_shared_library = '/usit/abel/u1/oleksanf/precimed/BGMG/src/build/lib/libbgmg.so';
bgmg_shared_library_header = '/usit/abel/u1/oleksanf/precimed/BGMG/srcbgmg_matlab.h';
logfile = 'BGMG_cpp_example.bgmglib.log';
BGMG_cpp.unload(); 
BGMG_cpp.load(bgmg_shared_library, bgmg_shared_library_header);
BGMG_cpp.init_log(logfile);
```
Now inspect logfile to see that there is "a new session" string in it. BGMG_cpp is a matlab class defined in BGMG_cpp.m file in the root of precimed/BGMG repository.


** Build on lisa.surfsara.nl

```
module purge
module load gcc/5.4.0
module load cmake/3.5.1

cd ~ && wget http://sourceforge.net/projects/boost/files/boost/1.49.0/boost_1_49_0.tar.gz && tar -xzvf boost_1_49_0.tar.gz && rm boost_1_49_0.tar.gz
./bootstrap.sh --with-libraries=program_options,filesystem,system,date_time i
./b2 cxxflags=-fPIC cflags=-fPIC variant=release threading=multi runtime-link=shared link=shared -j12
./b2 cxxflags=-fPIC cflags=-fPIC variant=release threading=multi runtime-link=static link=static -j12
mkdir ~/GitHub/precimed && cd ~/GitHub/precimed
git clone --recurse-submodules git@github.com:precimed/BGMG.git && cd BGMG/src
mkdir build && cd build
cmake .. -DBOOST_ROOT=/home/oleksand/boost_1_49_0
make -j12

```

** Run on lisa.surfsara.nl

```
module load matlab/2015b
export "LD_LIBRARY_PATH=/home/oleksand/GitHub/precimed/BGMG/src/build/lib:$LD_LIBRARY_PATH"
cd /home/oleksand/GitHub/precimed/BGMG

matlab -nodisplay -nosplash -nodesktop -r "
 trait1_file='/home/oleksand/GitHub/precimed/BGMG_reference/LDSR/LDSR_Data/PGC_SCZ_2014_noMHC.sumstats.gz';
 bim_file='/home/oleksand/GitHub/precimed/BGMG_reference/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim';
 frq_file='/home/oleksand/GitHub/precimed/BGMG_reference/LDSR/1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.@.frq';
 plink_ld_bin='/home/oleksand/GitHub/precimed/BGMG_reference/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin'; chr_labels = 1:22;
 out_file='/home/oleksand/GitHub/precimed/BGMG_reference/BGMG_result/PGC_SCZ_2014_noMHC.model=full.r2min=p05.randprune=n64p05.kmax=5000.run1.fit';
 bgmg_shared_library='src/build/lib/libbgmg.so'; bgmg_shared_library_header='src/bgmg_matlab.h';
 kmax=5000; max_causal_fraction=0.03; cache_tag_r2sum=0; SEED=123; randprune_r2=0.1; randprune_n=64; CI_ALPHA=0.05; r2min=0.0;
 DO_FIT_UGMG=1; FIT_FULL_MODEL=1; POWER_PLOT=1; POWER_PLOT_DOWNSCALE=100; QQ_PLOT=1; QQ_PLOT_DOWNSCALE=100; QQ_PLOT_BINS=1; QQ_PLOT_BINS_DOWNSCALE=50;
 UGMG_cpp_run_simple; "
```
The argument to matlab call has to be in one string; or, start matlab and enter all commands that set parameters, followed by `UGMG_cpp_run_simple` to trigger the actual run.

To update access: ``setfacl -R -m u:<user>:rx /home/oleksand/GitHub/precimed``.


**Building for MiXeR python wrapper on Abel**

```
module load cmake/3.7.1 python2/2.7.15.gnu
wget https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz 
tar -xzvf boost_1_69_0.tar.gz && cd boost_1_69_0
./bootstrap.sh --with-libraries=program_options,filesystem,system,date_time && ./b2 --clean && ./b2 --j12 -a
cmake .. -DBOOST_ROOT=/usit/abel/u1/oleksanf/boost_1_69_0 
make -j16 
```

I used ``python2/2.7.15.gnu`` to compile boost and mixer, but then run with ``python3/3.7.0.gnu``.
The toolchains are almost equal so it should make no difference when re-compile with ``python3/3.7.0.gnu``.
```
  1) binutils/2.26        2) gcc/6.1.0            3) openmpi.gnu/1.10.6   4) python2/2.7.15.gnu
  1) binutils/2.26        2) gcc/6.1.0            3) openmpi.gnu/1.10.6   4) openssl/1_1_1        5) python3/3.7.0.gnu
```

**Build on XSEDE Comet**

```
export MODULEPATH=$MODULEPATH:/share/apps/compute/modulefiles
module purge && module load gnutools/2.69 cmake/3.9.1 gnu/7.2.0  # gnutools/2.69 must go before gnu/7.2.0
cmake -DBOOST_ROOT=/home/oleksanf/boost_1_69_0 -DCMAKE_C_COMPILER="$(which gcc)" -DCMAKE_CXX_COMPILER="$(which g++)" ..
```

**Build on TSD**

```
# download CMake and Boost, build it from source
module load python3.gnu/3.7.3

cd /cluster/projects/p33/users/ofrei/no-backup/software/cmake-3.15.5
./bootstrap --prefix=/cluster/projects/p33/users/ofrei && make && make install

cd /cluster/projects/p33/users/ofrei/no-backup/software/boost_1_69_0
./bootstrap.sh --with-libraries=program_options,filesystem,system,date_time && ./b2 --clean && ./b2 --j12 -a
```
