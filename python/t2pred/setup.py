from distutils.core import setup, Extension
from Cython.Distutils import build_ext

setup(
        name='t2pred',
        version='1.0',
        description='Interface to tempo2 predictor library',
        author='M. Keith',

    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("t2pred", ["t2pred.pyx"],libraries=['tempo2pred'])]
)
