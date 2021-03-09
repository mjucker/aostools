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

## xarray integration

While my goal is to make all functions work with [xarray](xarray.pydata.org), this is still work in progress. For now, some of the functions exist for both numpy and xarray data, and if you have xarray Datasets or DataArrays, look for functions named `SomeFunctionNameXr()`, whereas the numpy equivalend would be `SomeFunctionName()`. For instance, the Eliassen-Palm flux calculations are done in `ComputeEPfluxDiv()` for numpy arrays and `ComputeEPfluxDivXr()` for xarray.DataArrays.

## How to cite

If you use any of the `aostools` functionality for your published work, please include a citation using either the generic DOI for all versions, [10.5281/zenodo.597598](https://doi.org/10.5281/zenodo.597598), or the DOI linking to the specific release, which you can find by visiting [the same link](https://doi.org/10.5281/zenodo.597598).

If you use the Eliassen-Palm flux calculations or plotting abilities, please also cite [Jucker ASL (2021), DOI 10.1002/asl.1020](https://doi.org/10.1002/asl.1020). If you use the wave activity flux calculations, please cite [Takaya & Nakamura GRL (1997), DOI 10.1029/97GL03094](https://doi.org/10.1029/97GL03094).

Finally, if you use any of the xarray capability, you might want to thank those developers by citing [Hoyer & Hamman JORS (2017), DOI 10.5334/jors.148](https://doi.org/10.5334/jors.148).
