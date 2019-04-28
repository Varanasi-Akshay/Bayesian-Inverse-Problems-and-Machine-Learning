from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from glob import glob
import numpy as np
import os

os.environ["CC"]="clang++"
os.environ["CXX"]="clang++"

exts  =  [Extension( "pyStructModel", glob('*.pyx') + glob('*.cpp'), include_dirs=[np.get_include()], language="c++", extra_compile_args=['-std=gnu++11']), ]

setup(
   name = "pyStructModel",
   ext_modules = cythonize( exts )
   )
