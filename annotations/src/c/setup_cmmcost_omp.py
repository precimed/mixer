# python setup_cmmcost_omp.py build_ext --inplace
# export OMP_NUM_THREADS=6
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy

srcfiles = ["cmmcost_omp.pyx", "_cmmcost_omp.c"]
includedirs = ["/mnt/seagate10/gsl/include", numpy.get_include()] # path to GSL include dir
libraries = ["gsl","gslcblas","m"]
libdirs = ["/mnt/seagate10/gsl/lib"]
compile_args = ["-O2", "-ffast-math", "-fopenmp"]
link_args = [] # "-flto"

ext = Extension("cmmcost_omp", sources=srcfiles, include_dirs=includedirs,
    library_dirs=libdirs, libraries=libraries, extra_link_args=link_args, 
    extra_compile_args=compile_args)

setup(ext_modules = cythonize([ext]))
