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

