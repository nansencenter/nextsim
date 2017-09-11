# compile bamg.so with:
# python setup.py build_ext --inplace

from distutils.core import setup, Extension
from Cython.Build import cythonize
import os


# ========================================================================
# source files
sources = [] # .cpp
pyxsrc  = [] # .pyx files
pyxsrc_ = [] # .cpp files from .pyx files
for src in os.listdir('.'):
   if len(src)>4:
      if src[-4:]=='.cpp':
         sources.append(src)
      elif src[-4:]=='.pyx':
         pyxsrc.append(src)
         pyxsrc_.append(src[:-4]+'.cpp')

# remove .cpp files made from 
for cppfil in pyxsrc_:
   if cppfil in sources:
      sources.remove(cppfil)

sources.extend(pyxsrc)
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
