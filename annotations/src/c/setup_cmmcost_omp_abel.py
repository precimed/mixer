# python3 setup_cmmcost_omp_abel.py build_ext --inplace
# export OMP_NUM_THREADS=16
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy

srcfiles = ["cmmcost_omp.pyx", "_cmmcost_omp.c"]
includedirs = ["/cluster/software/VERSIONS/gsl-2.2/include", numpy.get_include()] # path to GSL include dir
libraries = ["gsl","gslcblas","m"]
libdirs = ["/cluster/software/VERSIONS/gsl-2.2/lib"]
compile_args = ["-std=c99", "-O2", "-no-prec-div", "-qopenmp"]
link_args = ["-qopenmp"] # "-flto"

ext = Extension("cmmcost_omp", sources=srcfiles, include_dirs=includedirs,
    library_dirs=libdirs, libraries=libraries, extra_link_args=link_args, 
    extra_compile_args=compile_args)

setup(ext_modules = cythonize([ext]))
