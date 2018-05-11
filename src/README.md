Cost function for Univariate and Bivariate mixture model
Written as a plain C interface. Implemented in C++. Require Boost libraries.


To build on Windows:
- install MS VS 2015
- install pre-built boost libraries from https://sourceforge.net/projects/boost/files/boost-binaries/
- set env variables: 
```
setx BOOST_ROOT H:\local\boost_1_62_0
setx BOOST_LIBRARYDIR H:\local\boost_1_62_0\lib64-msvc-14.0
```
- ``cmake .. -G"Visual Studio 14 2015 Win64"``
- open bgmg.sln in Visual Studio, and build as usual
