from distutils.core import setup
from Cython.Build import cythonize

setup(
    name="testfuck",
    ext_modules = cythonize("TubeFunctions.pyx")
)