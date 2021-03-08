# aostools 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.597598.svg)](https://doi.org/10.5281/zenodo.597598) [![pypi](https://badge.fury.io/py/aostools.svg)](https://badge.fury.io/py/aostools)

Helper functions for scientific postprocessing and analysis of netCDF data. Also includes I/O routines to seamlessly work with [pv_atmos](https://github.com/mjucker/pv_atmos). 
This readme will be extended in the near future. Until then, each function has documentation, just try `help(myfunction)` for information.

## Installation

`aostools` is on PyPi. Thus, simply install the package running
```
pip install aostools
```
If you want to use map projections (with the `Projection()` function), you'd also want to install [Cartopy](https://scitools.org.uk/cartopy/docs/latest/). The best way to do this is via `conda`
```
conda install cartopy
```

