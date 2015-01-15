 #!/usr/bin/python
# Filename: tools.py
#
# Code by Martin Jucker, distributed under an MIT License
# Any publication benefitting from this piece of code should cite
# Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView.
# Journal of Open Research Software 2(1):e4, DOI: http://dx.doi.org/10.5334/jors.al
#
# This file provides helper functions that can be useful as pre-viz-processing of files and data

############################################################################################

## check if dimension already exists, and if so, create new dimension name
#
def MakeNewDimensionName(file,dimName):
    """Check if desired dimension name already exists.
       If so, add number until dimension name is unique

       INPUTS:
       file:    netCDF4 dataset
       dimName: desired dimension name
       OUPUT:
       newName: new unique dimension name
    """
    if dimName not in file.dimensions:
        newName = dimName
    else:
        newName = dimName
        i = 0
        while newName in file.dimensions:
            i += 1
            newName = dimName + str(i)
        return newName

## check if the given dimension already exists in a file. 
#
def CheckDimension(file,dimName,dimVal):
    """Check if desired dimension already exists in a file.
       If not, returns new dimension name.

       INPUTS:
       file:    netCDF4 dataset
       dimName: dimension name to be checked
       dimVal:  value of dimension to be checked
       OUTPUT:
       'exists':   dimension already exists and is identical
       'someName': dimension does not exist - returns unique new name
    """
    from numpy import allclose
    if dimName in file.dimensions:
        if dimName not in file.variables:
            raise Exception('Dimension '+dimName+' has no variable assiciated with it')
        dimExist = file.variables[dimName][:]
        if  allclose(dimExist,dimVal):
            return 'exists'
        else:
            return MakeNewDimensionName(file,dimName)
    else:
        return dimName

## create a generic netCDF file that can be read with pv_atmos
#
def WriteGenericNc(x, y, z, t, data, varName, outFileName='ncGeneric.nc'):
    """Write netCDF file that is compatible with pv_atmos.
        
        Takes one variable, and stores it inside a generic netCDF file.
        If the file already exists, adds the variable to the file.
        
        INPUTS:
        x,y,z,t     - spatial and time coordinates. If not all used, set to []
        data        - data to be written. Should be one variable.
        varName     - name of the variable data in the new netCDF file
        outFileName - name of (new) netCDF file. If exists already, variable is added.
        NOTE: if outFileName already exists, the dimensions in the file have to
        correspond to the dimensions given as inputs.
        """
    import netCDF4 as nc
    import os.path as op
    
    dims = ()
    if op.isfile(outFileName):
        mode = 'r+'
        text = 'Added to'
    else:
        mode = 'w'
        text = 'Created'
    outFile = nc.Dataset(outFileName, mode, format='NETCDF3_64BIT')
    # check if dimensions x,y,z already exist, and add them if not
    if len(x) > 0:
        checkDim = CheckDimension(outFile,'x',x)
        if checkDim == 'exists':
            checkDim = 'x'
        else:
            xD = outFile.createDimension(checkDim,len(x))
            xV = outFile.createVariable(checkDim,'f4', (checkDim,))
            xV[:] = x
            xV.setncattr('long_name','X axis')
            xV.setncattr('cartesian_axis','X')
        dims = (checkDim,) + dims
    if len(y) > 0:
        checkDim = CheckDimension(outFile,'y',y)
        if checkDim == 'exists':
            checkDim = 'y'
        else:
            yD = outFile.createDimension(checkDim,len(y))
            yV = outFile.createVariable(checkDim,'f4', (checkDim,))
            yV[:] = y
            yV.setncattr('long_name','Y axis')
            yV.setncattr('cartesian_axis','Y')
        dims = (checkDim,) + dims
    if len(z) > 0:
        checkDim = CheckDimension(outFile,'z',z)
        if checkDim == 'exists':
            checkDim = 'z'
        else:
            zD = outFile.createDimension(checkDim,len(z))
            zV = outFile.createVariable(checkDim,'f4', (checkDim,))
            zV[:] = z
            zV.setncattr('long_name','Z axis')
            zV.setncattr('cartesian_axis','Z')
        dims = (checkDim,) + dims
    if len(t) > 0:
        checkDim = CheckDimension(outFile,'time',t)
        if checkDim == 'exists':
            checkDim = 'time'
        else:
            tD = outFile.createDimension(checkDim,len(t))
            tV = outFile.createVariable(checkDim,'f4', (checkDim,))
            tV[:] = t
            tV.setncattr('long_name','time')
            tV.setncattr('cartesian_axis','T')
            tV.setncattr('units','time since 0000-00-00 00:00:00')
        dims = (checkDim,) + dims
    vOut = outFile.createVariable(varName,'f4',dims)
    vOut[:] = data
    outFile.close()
    print text+' file '+outFileName


