from Cython.Build import cythonize, build_ext
from numpy import get_include
from setuptools import setup, Extension, find_packages
import os
import numpy as np

extensions = [
    Extension("cythLeastR", ['pySPaRTAN/cythLeastR/cythLeastR.pyx', 'cythLeastR/ep21R.c', 'cythLeastR/epph.c'], include_dirs=['.', 'cythLeastR', get_include()]),
    Extension("cythKronPlus", ["pySPaRTAN/cythKronPlus/cythKronPlus.pyx"]),
    ]


setup(
    name='pySPaRTAN',
    ext_package=cythonize(extensions),
    zip_safe=False,
    packages=find_packages(),
    version="0.0.4",
    include_dirs=get_include(),
    cmdclass={'build_ext': build_ext}
)
