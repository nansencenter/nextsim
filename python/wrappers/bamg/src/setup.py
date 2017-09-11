# compile bamg.so with:
# python setup.py build_ext --inplace

from distutils.core import setup, Extension
from Cython.Build import cythonize
import os


# ========================================================================
# source files
sources = [] # .cpp, .pyx (need to be converted to .cpp)
for src in os.listdir('.'):
   if len(src)>4:
      if src[-4:]=='.cpp' or src[-4:]=='.pyx':
         sources.append(src)
# ========================================================================


# ========================================================================
# include dirs
nsdir    = os.getenv('NEXTSIMDIR')
incdirs  = []
incdirs.append('.')
incdirs.append('../include')
incdirs.append(nsdir+'/contrib/bamg/include')
# ========================================================================


# ========================================================================
# libraries used
libs     = []
libdirs  = []

# - bamg
libs   .append('bamg')
libdirs.append(nsdir+'/lib/')
# ========================================================================


# ========================================================================
# extra compilation arguments
xargs = ['-fopenmp'] # openMP
# ========================================================================


# ========================================================================
# compile
ext = Extension('lib/bamg',
                language='c++',
                sources=sources,
                include_dirs=incdirs,
                libraries=libs,
                library_dirs=libdirs,
                runtime_library_dirs=libdirs,
                extra_compile_args=xargs,
                extra_link_args=xargs
                )
setup(name='lib/bamg', ext_modules = cythonize(ext,build_dir='build'))
# ========================================================================
