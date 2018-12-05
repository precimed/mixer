# python setup_cmmcost.py build_ext --inplace
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

srcfiles = ["cmmcost.pyx", "_cmmcost.c"]
includedirs = ["/home/alexeas/.local/include"] # path to GSL include dir
libraries = ["gsl","gslcblas","m"]
libdirs = ["/home/alexeas/.local/lib"]
compile_args = ["-ffast-math"]

ext = Extension("cmmcost", sources=srcfiles, include_dirs=includedirs,
    library_dirs=libdirs, libraries=libraries, extra_compile_args=compile_args)

setup(ext_modules = cythonize([ext]))