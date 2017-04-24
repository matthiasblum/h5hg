#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension
from Cython.Build import cythonize

setup(
    name='h5hg',
    version='1.0',
    description='h5hg is a module to retrieve NGS data stored in the HDF5 format',
    author='Matthias Blum',
    author_email='mat.blum@gmail.com',
    url='https://github.com/matthiasblum/h5hg',
    zip_safe=False,
    install_requires=['numpy', 'scipy', 'h5py'],
    ext_modules=cythonize([
        Extension('h5hg', ['h5hg.py'])
    ])
)
