#from distutils.core import setup
import os
from setuptools import setup

def read(fname):
	return open(os.path.join(os.path.dirname(__file__),fname)).read()

setup(name='aostools',
      version='2.1.2',
      description='Helper functions for postprocessing and analysis of netCDF data',
      long_description=read('README.md'),
      author='Martin Jucker',
      author_email='coding@martinjucker.com',
      license='GPLv3',
      url='https://github.com/mjucker/aostools',
      package_dir={'aostools': ''},
      packages=['aostools'],
      install_requires=[
	'numpy',
	'netCDF4',
	],
      zip_safe=False
      )
