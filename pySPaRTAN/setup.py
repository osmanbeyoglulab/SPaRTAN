from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from numpy import get_include

import os
# rename exist cython extension files
cur_dir = os.path.dirname(os.path.abspath(__file__))

for filename in os.listdir(cur_dir):
    base_file, ext = os.path.splitext(filename)
    if base_file[-4:] == "_old":
        os.remove(os.path.join(cur_dir,filename) )
        
for filename in os.listdir(cur_dir):
    base_file, ext = os.path.splitext(filename)        
    if (ext == ".pyd" or ext == ".os"):                          
        os.rename(filename, base_file + "_old" + ext)

cur_file_ckron = os.path.join(cur_dir,'cythKronPlus','cythKronPlus.c')
old_file_ckron = os.path.join(cur_dir,'cythKronPlus','cythKronPlus_old.c')
cur_file_cleastR = os.path.join(cur_dir,'cythLeastR','cythLeastR.c')
old_file_cleastR = os.path.join(cur_dir,'cythLeastR','cythLeastR_old.c')

if os.path.exists(old_file_ckron):
    os.remove(old_file_ckron)
if os.path.exists(old_file_cleastR):
    os.remove(old_file_cleastR)
if os.path.exists(cur_file_ckron):
    os.rename(cur_file_ckron, old_file_ckron)
if os.path.exists(cur_file_cleastR):
    os.rename(cur_file_cleastR, old_file_cleastR)

extensions = [
    Extension("cythLeastR", ['cythLeastR/cythLeastR.pyx', 'cythLeastR/ep21R.c', 'cythLeastR/epph.c'], include_dirs=['.', 'cythLeastR', get_include()]),
    Extension("cythKronPlus", ["cythKronPlus/cythKronPlus.pyx"]),
    ]

setup(
    ext_modules=cythonize(extensions),
)


