#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(
    name='h5hg',
    version='1.0',
    description='h5hg is a module to retrieve NGS data stored in the HDF5 format',
    author='Matthias Blum',
    author_email='mat.blum@gmail.com',
    url='https://github.com/matthiasblum/h5hg',
    py_modules=['h5hg'],
    zip_safe=False,
    install_requires=['numpy', 'scipy', 'h5py']
)
