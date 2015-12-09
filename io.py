 #!/usr/bin/python
# Filename: io.py
#
# Code by Martin Jucker, distributed under an MIT License
# Any publication benefitting from this piece of code should cite
# Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView.
# Journal of Open Research Software 2(1):e4, DOI: http://dx.doi.org/10.5334/jors.al
#
# This file provides helper functions that can be useful as pre-viz-processing of files and data

############################################################################################


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

## create a generic netCDF file that can be read with pv_atmos
#
def WriteGenericNc(x, y, z, t, data, varName, outFileName='ncGeneric.nc',dimNames=['x','y','z','time'],unlimTime=True,timeUnits='days since 0001-01-01 00:00:00'):
    """Write netCDF file that is compatible with pv_atmos.

        Note that due to differences bewteen python and netCDF dimensions,
        data.shape = (len(t),len(z),len(y),len(x))
        
        Takes one variable, and stores it inside a generic netCDF file.
        If the file already exists, adds the variable to the file.
        
        INPUTS:
        x,y,z,t     - spatial and time coordinates. If not all used, set to []
        data        - data to be written. Should be one variable.
        varName     - name of the variable data in the new netCDF file
        outFileName - name of (new) netCDF file. If exists already, variable is added.
        dimNames    - names of the output dimensions
        unlimTime   - whether time dimension should be unlimited/appendable
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
        checkDim = CheckDimension(outFile,dimNames[0],x)
        if checkDim == 'exists':
            checkDim = dimNames[0]
        else:
            xD = outFile.createDimension(checkDim,len(x))
            xV = outFile.createVariable(checkDim,'f4', (checkDim,))
            xV[:] = x
            xV.setncattr('long_name','X axis')
            xV.setncattr('cartesian_axis','X')
        dims = (checkDim,) + dims
    if len(y) > 0:
        checkDim = CheckDimension(outFile,dimNames[1],y)
        if checkDim == 'exists':
            checkDim = dimNames[1]
        else:
            yD = outFile.createDimension(checkDim,len(y))
            yV = outFile.createVariable(checkDim,'f4', (checkDim,))
            yV[:] = y
            yV.setncattr('long_name','Y axis')
            yV.setncattr('cartesian_axis','Y')
        dims = (checkDim,) + dims
    if len(z) > 0:
        checkDim = CheckDimension(outFile,dimNames[2],z)
        if checkDim == 'exists':
            checkDim = dimNames[2]
        else:
            zD = outFile.createDimension(checkDim,len(z))
            zV = outFile.createVariable(checkDim,'f4', (checkDim,))
            zV[:] = z
            zV.setncattr('long_name','Z axis')
            zV.setncattr('cartesian_axis','Z')
        dims = (checkDim,) + dims
    if len(t) > 0:
        checkDim = CheckDimension(outFile,dimNames[3],t)
        if checkDim == 'exists':
            checkDim = dimNames[3]
        else:
            if unlimTime:
                tlen = 0
            else:
                tlen = len(t)
            tD = outFile.createDimension(checkDim,tlen)
            tV = outFile.createVariable(checkDim,'f4', (checkDim,))
            tV[:] = t
            tV.setncattr('long_name','time')
            tV.setncattr('cartesian_axis','T')
            tV.setncattr('units',timeUnits)
        dims = (checkDim,) + dims
    vOut = outFile.createVariable(varName,'f4',dims)
    vOut[:] = data
    outFile.close()
    print text+' file '+outFileName
    print 'Used dimensions '+str(dims)
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
        print 'Done, written file '+outFileName
    else:
        print 'Done, appended variable to '+outFileName

## Read a file, and show the variables contained in it
#
def ReadFile(fileName, show='short', mode='r'):
    """Reads netCDF4 file fileName and shows available variables.
        
        INPUTS:
        fileName: full path to file
        show:     how much detail to show; 'min' for nothing, 'short' for one line, 'full'
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
    if show != 'min':
        print 'All variables:',file.variables.keys()
    if show == 'full':
        for v in file.variables.keys():
            if v not in file.dimensions:
                print file.variables[v]
    return file

