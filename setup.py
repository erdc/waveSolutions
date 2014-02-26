#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy,mpi4py

setup(name='waveSolutions',
      version='0.0.1',
      description='Wave modules for analytical  water waves.',
      author='Matt Malej',
      author_email='matt.malej@erdc.dren.mil',
      url='https://github.com/erdc-cm/waveSolutions',
      packages=['randomTests','randomTests.testWaves'],
      package_dir={'testWaves':'randomTests/testWaves'},
      package_data={'randomTests':['testWaves/*.py']},
      scripts = ['runWaves.py'],
      requires=['numpy']
      )
