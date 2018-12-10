from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import shutil
import numpy as np
ext_modules = [Extension("gem_qa", ["gem_qa.pyx"], include_path = [np.get_include(), './'], include_dirs = [np.get_include(), './'], 
                         language = 'c++', extra_compile_args=["-std=c++11"], extra_link_args=["-std=c++11"]),\
               Extension("CSeqDict", ["CSeqDict.pyx"], include_path = ['./'], include_dirs = [], language = 'c++', extra_compile_args=["-std=c++11"],
                         extra_link_args=["-std=c++11"]),
               Extension("Delta", ["Delta.pyx"])]
setup(cmdclass={'build_ext': build_ext}, ext_modules=ext_modules)
shutil.copy('build/lib.linux-x86_64-2.7/gem_qa.so', 'gem_qa.so')
shutil.copy('build/lib.linux-x86_64-2.7/CSeqDict.so', 'CSeqDict.so')
shutil.copy('build/lib.linux-x86_64-2.7/Delta.so', 'Delta.so')