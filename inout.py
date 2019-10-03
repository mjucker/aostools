 #!/usr/bin/python
# Filename: inout.py
#
# Code by Martin Jucker, distributed under an GPLv3 License
#
# This file provides helper functions that can be useful as pre-viz-processing of files and data

############################################################################################
from __future__ import print_function
## find compression parameters as in ncpdq
def DefCompress(x,varName=None):
    """Produce encoding dictionary for to_netcdf(encoding=encodeDict).

        INPUTS:
            x       : either a xarray.DataArray or xarray.Dataset
            varName : only encode that variable. If None, encode all variables
        OUTPUTS:
            encodeDict : dictionary containing compressing information for each
                          variable. To be used as x.to_netcdf(encoding=encodeDict)
    """
    import numpy as np
    # check whether x is a Dataset or DataArray
    try:
        keys = x.variables.keys()
        is_dataset = True
    except:
        is_dataset = False
    # make a list of variables
    if varName is None:
        vars = []
        if is_dataset:
            for var in x.variables.keys():
                vars.append(var)#.encode("utf-8"))
        else:
            vars.append(x.name)#.encode("utf-8"))
    else:
        if isinstance(varName,list):
            vars = varName
        else:
            vars = [varName]
    # now loop over all variables
    encodeDict = {}
    bytes = 16
    fillVal = -2**(bytes-1)
    for var in vars:
        if is_dataset:
            filtr = np.isfinite(x[var])
            dataMin = x[var].where(filtr).min().values
            dataMax = x[var].where(filtr).max().values
        else:
            filtr = np.isfinite(x)
            dataMin = x.where(filtr).min().values
            dataMax = x.where(filtr).max().values
        scale_factor=(dataMax - dataMin) / (2**bytes - 2)
        add_offset = (dataMax + dataMin) / 2
        encodeDict[var] = {
            'dtype':'short',
            'scale_factor':scale_factor,
            'add_offset': add_offset,
            '_FillValue': fillVal}
    return encodeDict

## copy all dimensions from one file to the other
def CopyDims(inFile, outFile, onlyDims=[]):
    """Copy all dimensions from one file to the other.

        INPUTS:
        inFile  - netCDF4 Dataset of file containing the dimensions
        outFile - netCDF4 Dataset of file to which dimensions will be copied
        onlyDims- copy only specified dimensions
        OUTPUTS:
        outFile - netCDF4 Dataset containing the same dimensions as inFile
        """
    for dim in inFile.dimensions:
        copyDim = True
        if len(onlyDims) > 0 and dim not in onlyDims:
            copyDim = False
        if copyDim:
            outDim = outFile.createDimension(dim,len(inFile.dimensions[dim]))
            inVar = inFile.variables[dim]
            outVar = outFile.createVariable(dim, str(inFile.variables[dim].dtype),(dim,))
            outVar = CopyAtts(inVar, outVar)
            outVar[:] = inVar[:]
    return outFile
#
## copy all attributes of one variable to another
def CopyAtts(inVar, outVar):
    """Copy all attributes of one variable to the other

        INPUTS:
        inVar  - netCDF4 variable containing attributes to copy
        outVar - netCDF4 variable to which attributes will be copied
        OUTPUTS:
        outVar - netCDF4 variable with all copied attributes
        """
    for att in inVar.ncattrs():
        outVar.setncattr(att,inVar.getncattr(att))
    return outVar
#

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
            raise Exception('Dimension '+dimName+' has no variable associated with it')
        dimExist = file.variables[dimName][:]
        if  allclose(dimExist,dimVal):
            return 'exists'
        else:
            return MakeNewDimensionName(file,dimName)
    else:
        return dimName

def MakeTimeDim(timeVar,timeUnits,calendar=None):
    timeVar.setncattr('long_name','time')
    timeVar.setncattr('cartesian_axis','T')
    timeVar.setncattr('units',timeUnits)
    if calendar:
        timeVar.setncattr('calendar',calendar)
    return timeVar

## create a generic netCDF file that can be read with pv_atmos
#
def WriteGenericNc(dimensions, data, varName, outFileName='ncGeneric.nc',dimNames=None,timeDim=None,unlimTime=True,timeUnits='days since 0001-01-01 00:00:00',calendar=None,group=None,zlib=False,least_digit=None):
    """Write netCDF file that is compatible with pv_atmos.

        Note that due to differences bewteen python and netCDF dimensions,
        the dimensions are in reverse order. For instance,
        dimensions = (x,y,z,t)
        data.shape = (len(t),len(z),len(y),len(x))

        Takes one variable, and stores it inside a generic netCDF file.
        If the file already exists, adds the variable to the file.

        INPUTS:
        dimenions   - list of all spatial and time coordinates.
        data        - data to be written. Should be one variable.
        varName     - name of the variable data in the new netCDF file
        outFileName - name of (new) netCDF file. If exists already, variable is added.
        dimNames    - names of the output dimensions
        timeDim     - index of the time variable
        unlimTime   - whether time dimension should be unlimited/appendable
        timeUnits   - units of time dimension
        calendar    - calendar of time dimension
        group       - if not None, must be string with name of groupe to write to
        zlib        - use zlib compression
        least_digit - further compress using least_significant_digit
        NOTE: if outFileName already exists, the dimensions in the file have to
        correspond to the dimensions given as inputs.
        """
    import netCDF4 as nc
    import os.path as op

    letters = ['a','b','c','d','e','f','g','h','i','k','l','m','n']

    dims = ()
    if op.isfile(outFileName):
        mode = 'r+'
        text = 'Added to'
    else:
        mode = 'w'
        text = 'Created'
    # outFile = nc.Dataset(outFileName, mode, format='NETCDF3_64BIT')
    outFile = nc.Dataset(outFileName, mode, format='NETCDF4')
    # check if group already exists, if not, create
    if group is not None:
        if group not in outFile.groups.keys():
            outFile.createGroup(group)
        group = '/'+group+'/'
    else:
        group = ''
    # check if dimensions already exist, and add them if not
    for d in range(len(dimensions)):
        dim = dimensions[d]
        try:
            name = dimNames[d]
        except:
            name = letters[d]
        checkDim = CheckDimension(outFile,name,dim)
        if checkDim == 'exists':
            checkDim = name
        else:
            if timeDim == d and unlimTime:
                dlen = 0
            else:
                dlen = len(dim)
            xD = outFile.createDimension(checkDim,dlen)
            xV = outFile.createVariable(checkDim,'f4', (checkDim,))
            xV[:] = dim
            if timeDim == d:
                xV = MakeTimeDim(xV,timeUnits,calendar)
        dims = (checkDim,) + dims
    # add the variable
    if least_digit is None:
        vOut = outFile.createVariable(group+varName,'f4',dims,zlib=zlib)
    else:
        vOut = outFile.createVariable(group+varName,'f4',dims,zlib=zlib,least_significant_digit=least_digit)
    vOut[:] = data
    outFile.close()
    print( text+' file '+outFileName )
    print( 'Used dimensions '+str(dims) )
#
# write a new file, with the same structure as an existing file
#  or add data to existing file
def WriteFileLike(data, varName, outFileName, dimNames, inFileName):
    """Write a new file with data, but dimensions from an existing file.
        Appends the new variable if outFileName==inFileName.

        New file will have one field, data, with name varName, and the
        dimensions given in dimNames. Also takes care of edge dimensions.

        INPUTS:
        data        - the variable to be put into the new file
        varName     - the name of the variable in the file
        outFileName - the name of the new file
        dimNames    - names of the dimensions to be copied from
        the existing file
        inFileName  - name of the existing file to copy dimensions
        if inFileName == outFileName, append variable
        """
    import netCDF4 as nc
    if outFileName != inFileName:
        inFile  = nc.Dataset(inFileName ,'r')
        outFile = nc.Dataset(outFileName,'w',format=inFile.file_format)
        dimNamesTmp = dimNames[:]
        for dim in dimNames:
            try:
                edgeDim = inFile.variables[dim].getncattr('edges')
                dimNamesTmp.append(edgeDim)
            except:
                pass
        outFile = CopyDims(inFile, outFile, dimNamesTmp)
        inFile.close()
    else:
        outFile = nc.Dataset(inFileName ,'a')
    outVar = outFile.createVariable(varName, 'f8', dimNames)
    outVar[:] = data
    outFile.close()
    if outFileName != inFileName:
        print( 'Done, written file '+outFileName )
    else:
        print( 'Done, appended variable to '+outFileName )

## Read a file, and show the variables contained in it
#
def ReadFile(fileName, show='silent', mode='r'):
    """Reads netCDF4 file fileName and shows available variables.

        INPUTS:
        fileName: full path to file
        show:     how much detail to show; 'silent' for nothing, 'short' for one line, 'full'
                    for full dimensionality and attributes
        mode:     'r' for read only, 'w' for write, 'r+' for read-write
        OUTPUTS:
        file:     netCDF4 Dataset
    """
    import os
    import netCDF4 as nc
    if not os.path.isfile(fileName):
        raise IOError('File '+fileName+' does not exist')
    file=nc.Dataset(fileName,mode)
    if show != 'silent':
        print( 'All variables:',file.variables.keys() )
        for g in file.groups.keys():
            print( '              ','group '+g+':' )
            print( '              ',file.groups[g].variables.keys() )
    if show == 'full':
        for v in file.variables.keys():
            if v not in file.dimensions:
                print( file.variables[v] )
        for g in file.groups.keys():
            for v in file.groups[g].variables.keys():
                if v not in file.dimensions:
                    print( g+'/'+file.variables[v] )
    return file


#######################################################
def CombineFiles(fileList, fileName, combineAxis=0, combineDim='time', meanAxis=None, meanDim=None, fileFormat='NETCDF3_64BIT', compress=False):
    """Combine a number of files into one along a given dimension

       The option to also average out one dimension does not work yet.
       For ParaView, the fileFormat should probably be NETCDF3_64BIT, but
       for compression, it has to be NETCDF4 or NETCDF4_CLASSIC.
    """
    import numpy as np
    import netCDF4 as nc


    # create output file
    outFile = nc.Dataset(fileName,'w',format=fileFormat)

    # copy all dimensions from first file in list to new file
    file = nc.Dataset(fileList[0],'r')
    for dim in file.dimensions:
        if dim == combineDim:
            dlen = 0 # unlimited
        else:
            dlen = len(file.dimensions[dim])
        outDim = outFile.createDimension(dim,dlen)
        inVar  = file.variables[dim]
        outVar = outFile.createVariable(dim,str(file.variables[dim].dtype),(dim,))
        outVar = CopyAtts(inVar,outVar)
        outVar[:] = inVar[:]

    # now add all variables
    outVars = {}
    outVars[combineDim] = outFile.variables[combineDim]
    for var in file.variables:
        if not var in outFile.variables:
            dims = file.variables[var].dimensions
            # remove the averaging dimension
            if meanAxis is not None and var is not combineDim:
                dims = tuple(d for d in dims if d != meanDim)
            outVars[var] = outFile.createVariable(var,'f4',dims,zlib=compress,complevel=9)
            inVar = file.variables[var]
            for att in inVar.ncattrs():
                if att != 'scale_factor' and att != 'add_offset':
                    outVars[var].setncattr(att,inVar.getncattr(att))
    # get all variables from first file
    tmpVars = {}
    for var in [combineDim]+outVars.keys():
        tmpVar = file.variables[var][:]
        if meanAxis is not None:
            tmpVar = tmpVar.mean(axis=meanAxis)
        tmpVars[var] = tmpVar[:]
    file.close()
    print( 'done with file',fileList[0] )

    # finally, add all other files
    for f in fileList[1:]:
        file = nc.Dataset(f,'r')
        for var in outVars.keys():
            tmpVar = file.variables[var][:]
            if meanAxis is not None:
                tmpVar = tmpVar.mean(axis=meanAxis)
            if var == combineDim:
                tmpVars[var] = np.concatenate((tmpVars[var],tmpVar))
            else:
                tmpVars[var] = np.concatenate((tmpVars[var],tmpVar),axis=combineAxis)
        file.close()
        print( 'done with file',f )
    #
    # and put all of it into the new file
    outFile.variables[combineDim][:] = tmpVars[combineDim][:]
    for var in outVars.keys():
        outVars[var][:] = tmpVars[var][:]
    outFile.close()
