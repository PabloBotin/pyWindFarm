from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    name='Gauss Seidel Iteration',
    ext_modules=cythonize("gauss_seidel.pyx"),
    include_dirs=[numpy.get_include()],
    extra_compile_args=["-O3"],
)