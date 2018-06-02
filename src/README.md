Cost function for Univariate and Bivariate mixture model
Written as a plain C interface. Implemented in C++. Require Boost libraries.


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

ldd /cluster/software/VERSIONS/matlab/R2017a/bin/glnxa64/MATLAB
export "LD_LIBRARY_PATH=/cluster/software/VERSIONS/matlab/R2017a/sys/os/glnxa64:$LD_LIBRARY_PATH"
ldd libbgmg.so
Errors like ``cluster/software/VERSIONS/matlab/R2017a/sys/os/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found`` indicate that you've built libbgmg.so with too new version of gcc.
