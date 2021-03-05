#from distutils.core import setup
import os
from setuptools import setup

def read(fname):
	return open(os.path.join(os.path.dirname(__file__),fname)).read()

setup(name='aostools',
      version='2.3',
      description='Helper functions for scientific postprocessing and analysis of netCDF data',
      long_description=read('readme_pypi'),
      keywords='atmospheric oceanic science netcdf analysis tools',
      author='Martin Jucker',
      author_email='coding@martinjucker.com',
      license='GPLv3',
      url='https://github.com/mjucker/aostools',
      package_dir={'aostools': ''},
      packages=['aostools'],
      install_requires=[
	'numpy',
	'netCDF4',
        'xarray>=0.11.0',
	'scipy',
	'cartopy',
	],
      zip_safe=False
      )
