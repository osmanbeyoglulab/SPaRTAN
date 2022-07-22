from Cython.Build import cythonize
from numpy import get_include
from setuptools import setup, Extension
import os

extensions = [
    Extension("cythLeastR", ['pySPaRTAN/cythLeastR/cythLeastR.pyx', 'cythLeastR/ep21R.c', 'cythLeastR/epph.c'], include_dirs=['.', 'cythLeastR', get_include()]),
    Extension("cythKronPlus", ["pySPaRTAN/cythKronPlus/cythKronPlus.pyx"]),
    ]


setup(
    name='pySPaRTAN',
    ext_package=cythonize(extensions),
    zip_safe=False
)
