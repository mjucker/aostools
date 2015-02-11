#!/usr/bin/python
# Filename: climate.py
#
# Code by Martin Jucker, distributed under an MIT License
# Any publication benefitting from this piece of code should cite
# Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView.
# Journal of Open Research Software 2(1):e4, DOI: http://dx.doi.org/10.5334/jors.al
#
# This file provides helper functions that can be useful as pre-viz-processing of files and data

############################################################################################
#
# compute climatologies
def ComputeClimate(file, climatType, wkdir='/', timeDim='time'):
    """Compute climatologies from netCDF files.
        
        ComputeClimate(file,climatType,wkdir='/',timeDim='time')
        
        Inputs:
        file        file name, relative path from wkdir
        climatType  'daily', 'monthly', 'DJF', 'JJA', 'annual', or any
                    combination of months according to two-letter code
                    Ja Fe Ma Ap My Jn Jl Au Se Oc No De
        wkdir       working directory, in with 'file' must be, and to which the output
                    is written
        timeDim     name of the time dimension in the netcdf file
        Outputs:
        writes outputfile with name depending on input file name and climatType
    """

    # need to read netCDF and of course do some math

    import netCDF4 as nc
    import numpy as np

    if climatType == 'DJF':
        climType = 'DeJaFe'
    elif climatType == 'JJA':
        climType = 'JuJlAu'
    elif climatType == 'annual':
        climType = 'JaFeMaApMyJnJlAuSeOcNoDe'
    else:
        climType = climatType
    monthList=['Ja','Fe','Ma','Ap','My','Jn','Jl','Au','Se','Oc','No','De']


    # Now, let's go and read that file. Also, find out how many time steps we have

    ncFile = nc.Dataset(wkdir+file,'r')

    numTimeSteps = len(ncFile.dimensions[timeDim])


    # Now, find out about the units and thus the time step interval in units of days

    timeVar = ncFile.variables[timeDim]
    if 'seconds' in timeVar.units:
        timeStepFact = 1./86400
        timeUnits = 'seconds'
    elif 'days' in timeVar.units:
        timeStepFact = 1
        timeUnits = 'days'
    elif 'months' in timeVar.units:
        timeStepFact = 30
        timeUnits = 'months'
    timeStep = np.diff(timeVar).mean()*timeStepFact
    print 'The time dimension is in units of',timeUnits,', with a time step of',timeStep,'days'


    # Now, make an educated guess about the length of the year in units of days.
    # At the moment, we can't deal with leap years.

    totalNumDays = ( timeVar[-1] - timeVar[0] ) * timeStepFact + timeStep
    if totalNumDays % 360 == 0:
        daysPerYear = 360
        daysPerMonths = 30*np.ones(12,)
    elif totalNumDays % 365 == 0:
        daysPerYear = 365
        daysPerMonths = [31,28,31,30,31,30,31,31,30,31,30,31]
    else:
        raise Exception("Cannot deal with years that are not 360 or 365 days. This data has "+str(totalNumDays)+" days.")
    stepsPerYear  = daysPerYear/timeStep
    totalNumYears = totalNumDays/daysPerYear
    print 'The simulation is',totalNumDays,'days long, which I assume means',totalNumYears,'years of',daysPerYear,'days.'


    # Now, need to know about the type of climatology we want.
    # 
    # climStep gives the number of time steps to include per climatology interval
    # 
    # climInts gives the indicies for each climatology interval

    # In[ ]:

    if   climType == 'daily':
        climStep = max(1,1/timeStep)
        climInts = climStep*np.ones(daysPerYear/max(1,timeStep),)
    elif climType == 'monthly':
        climInts = daysPerMonths/timeStep
    else:
        climInts = []
        monthEnds=np.cumsum(daysPerMonths)
        for m in range(len(monthList)):
            indx = (monthEnds[m]-monthEnds[m-1])%daysPerYear/timeStep
            if monthList[m] in climType:
                climInts.append(indx)
            else:
                climInts.append(-indx)


    # climInts is the important output here, as it knows for each year, which time steps will have to be averaged over
    # 
    # This is then put into includeSteps, which indicates which time steps one has to average over, considering the entire data set

    if climType == 'daily' or climType == 'monthly':
        includeSteps = np.zeros((stepsPerYear,len(climInts)))
        indx=0;m=0
        for n in climInts:
            includeSteps[indx:indx+n,m]=1
            indx = indx+n
            m += 1
    else:
        includeSteps = np.zeros(stepsPerYear,)
        indx=0
        for n in climInts:
            if n > 0:
                includeSteps[indx:indx+n]=1
            indx = indx+abs(n)
    includeSteps = np.array(list(includeSteps)*totalNumYears)


    # Create the output file, including dimensions.
    # 
    # We exclude time for seasonal climatologies, but need time for daily and monthly.

    outFileName = wkdir + file[0:-3] + '_' + climatType + '.nc'
    try:
        os.remove(outFileName)
    except:
        pass
    outFile = nc.Dataset(outFileName,'w',format=ncFile.file_format)
    for dim in ncFile.dimensions:
        if  dim != timeDim:
            outDim = outFile.createDimension(dim,len(ncFile.dimensions[dim]))
            inVar  = ncFile.variables[dim]
            outVar = outFile.createVariable(dim,str(ncFile.variables[dim].dtype),(dim,))
            outVar[:] = inVar[:]
            for att in inVar.ncattrs():
                outVar.setncattr(att,inVar.getncattr(att))
        elif climType == 'daily' or climType == 'monthly':
            nTime = np.shape(includeSteps)[1]
            if climType == 'daily':
                units = 'days'
                dTime = (np.arange(nTime) + 1)*max(1,timeStep)
            else:
                units = 'months'
                dTime = np.arange(nTime) + 1
            outDim = outFile.createDimension(dim,nTime)
            timeValue = dTime
            outVar = outFile.createVariable(dim,str(ncFile.variables[dim].dtype),(dim,))
            outVar[:] = timeValue
            outVar.setncattr('long_name','climatological ' + units[:-1] + ' of year')
            outVar.setncattr('units',units + ' since 0001-01-01 00:00:00')
            outVar.setncattr('cartesian_axis','T')
            outVar.setncattr('bounds','time_bounds')


    # Finally, perform the averaging and write into new file
    # 
    # Here, we need to be very careful in the event of packaged data: netCDF4 knows about packaging when reading data, but we need to use scale_factor and add_offset to package the data back when writing the new file.

    print 'Averaging variables:'
    for var in ncFile.variables:
        varShape = np.shape(ncFile.variables[var])
        if varShape[0] == numTimeSteps and len(varShape) >= 2:
            print '                     ',var
            tmpVar = ncFile.variables[var][:]
            if climType != 'daily' and climType != 'monthly':
                outVar = outFile.createVariable(var,str(ncFile.variables[var].dtype),ncFile.variables[var].dimensions[1:])
                tmpAvg = tmpVar[includeSteps>0,:].mean(axis=0)
            else:
                outVar = outFile.createVariable(var,str(ncFile.variables[var].dtype),ncFile.variables[var].dimensions    )
                avgShape = []
                avgShape.append(nTime)
                for t in range(len(np.shape(outVar))-1):
                    avgShape.append(np.shape(outVar)[t+1])
                tmpAvg = np.zeros(avgShape)
                for t in range(nTime):
                    tmpAvg[t,:] = tmpVar[includeSteps[:,t]>0,:].mean(axis=0)
            #package average
            if 'add_offset' in ncFile.variables[var].ncattrs():
                tmpAvg = tmpAvg - ncFile.variables[var].getncattr('add_offset')
            if 'scale_factor' in ncFile.variables[var].ncattrs():
                tmpAvg = tmpAvg/ncFile.variables[var].getncattr('scale_factor')
            #put the packaged average into the output variable
                outVar[:] = tmpAvg.astype(np.int16)
            else:
                outVar[:] = tmpAvg
            inVar = ncFile.variables[var]
            for att in inVar.ncattrs():
                outVar.setncattr(att,inVar.getncattr(att))

    ncFile.close()
    outFile.close()
    print 'DONE, wrote file',outFileName

##############################################################################################
def ComputeSaturationMixingRatio(T, p):
    """Computes the saturation water vapor mixing ratio according to Clausius-Clapeyron
        
        INPUTS:
            T - temperature in Kelvin
            p - pressure in hPa/mbar
        OUTPUTS:
            qsat - saturation water mixing ratio [kg/kg]
    """
    import numpy as np
    #some constants we need
    Rd = 287.04
    Rv = 461.5
    ES0 = 610.78
    HLV = 2.5e6
    Tf = 273.16
    
    # pressure is assumed in hPa: convert to Pa
    p = p*100
    # compute saturation pressure
    esat = ES0*np.exp(HLV*(1./Tf - 1./T)/Rv)
    # now, esat can have an arbitrary size
    #  we know that one dimension is the same for p and esat
    #  we also assume that all other dimensions have a different length
    I = find(np.array(np.shape(esat)) == len(p))
    # check that this dimension is unique
    if len(I) != 1:
        raise Exception("THE DIMENSIONS ARE NOT WELL DEFINED: CANNOT ASSIGN VERTICAL DIRECTION")
    # now construct an array containing p with the shape of esat
    sh = np.ones(len(np.shape(esat)))
    sh = sh/sh
    sh[I[0]] = len(p)
    pfull = np.zeros_like(p)
    pfull[:] = p
    for i in range(len(sh)):
        if sh[i] == 1:
            pfull = np.expand_dims(pfull,axis=i)
    # finally, compute saturation mixing ratio from pressure
    qsat = Rd/Rv*esat/(pfull-esat)
    return qsat


##############################################################################################
def ComputeRelativeHumidity(inFile, outFile='none', temp='temp', sphum='sphum', pfull='pfull', wkdir='/'):
    """Computes relative humidity from temperature and specific humidity.
        
        File inFile is assumed to contain both temperature and specific humidity.
        Relative humidity is either output of the function, or written to the file outFile.
        
        Inputs:
            inFile    Name of the file (path relative from wkdir)
                        containing temperature and moisture
            outFile   Name of the output file containing specific humidity.
                        No output file is created if outFile='none'
            temp      Name of the temperature variable inside inFile
            sphum     Name of specific humidity variable inside inFile
            pfull     Name of full level pressure [hPa] inside inFile
            wkdir     Path from which inFile and outFile is given relatively
    """

    import netCDF4 as nc
    import numpy as np
    # relative humidity is then q/qsat*100[->%]
    
    # read input file
    inFile = nc.Dataset(inFile, 'r')
    t = inFile.variables[temp]
    q = inFile.variables[sphum]
    p = inFile.variables[pfull]

    # compute saturation mixing ratio
    qsat = ComputeSaturationMixingRatio(t, p)

    #write output file
    if outFile is not 'none':
        outFile = nc.Dataset(inFile[0:-3]+'_out.nc','w')
        outFile = CopyDims(ncFile, outFile)
        outVar = outFile.createVariable('rh', 'f4', ncFile.variables[temp].dimensions)
        outVar[:] = q/qsat*1.e2
    return q/qsat*1.e2


##############################################################################################
def ComputePsi(inFileName, outFileName='same', temp='temp', vcomp='vcomp', lat='lat', pfull='pfull', time='time', p0=1e3, wkdir='./'):
    """Computes the residual stream function \Psi* (as a function of time).
        
        INPUTS:
            inFileName  - filename of input file, relative to wkdir
            outFileName - filename of output file, 'none', or 'same'
            temp        - name of temperature field in inFile
            vcomp       - name of meridional velocity field in inFile
            lat         - name of latitude in inFile
            pfull       - name of pressure in inFile [hPa]
            time        - name of time field in inFile. Only needed if outFile used
            p0          - pressure basis to compute potential temperature [hPa]
            wkdir       - root directory, from with the paths inFile and outFile are taken
        OUTPUTS:
            psi         - residual stream function, as a function of time
    """
    import netCDF4 as nc
    import numpy as np

    # some constants
    kappa = 2./7
    a0    = 6371000
    g     = 9.81

    # read input file
    print 'Reading data'
    update_progress(0)
    if outFileName == 'same':
        mode = 'a'
    else:
        mode = 'r'
    inFile = nc.Dataset(wkdir + inFileName, mode)
    t = inFile.variables[temp][:]
    update_progress(.25)
    v = inFile.variables[vcomp][:]
    update_progress(.50)
    l = inFile.variables[lat][:]
    update_progress(.75)
    p = inFile.variables[pfull][:]
    update_progress(1)

    ## compute psi
    print 'Computing psi'
    update_progress(0)

    v_bar,t_bar = compute_intv(v,t,p,p0) # t_bar = bar(v'Th'/(dTh_bar/dp))

    # Eulerian streamfunction
    dp  = np.gradient(p)[np.newaxis,:,np.newaxis]*100.
    psi = np.cumsum(v_bar*dp,axis=1)
    v_bar=v=t=[]
    update_progress(.50)


    ## compute psi* = psi - bar(v'Th'/(dTh_bar/dp))
    psis = psi - t_bar
    t_bar = []
    psi = 2*np.pi*a0/g*psi *np.cos(l[np.newaxis,np.newaxis,:]*np.pi/180.)
    psis= 2*np.pi*a0/g*psis*np.cos(l[np.newaxis,np.newaxis,:]*np.pi/180.)
    update_progress(1)

    ## write outputfile
    if outFileName is not 'none':
        print 'Writing psi* into file'
        if outFileName is not 'same':
            outFile = nc.Dataset(wkdir + outFileName,'w')
            outFile  = CopyDims(inFile, outFile, [time,pfull,lat])
        else:
            outFile = inFile
        outVar = outFile.createVariable('psi_star', 'f4', (time,pfull,lat,))
        outVar[:] = psis
        outFile.close()
        print 'Done writing file '+outFileName
        if outFileName is not 'same':
            inFile.close()
    return psi,psis


## helper functions
def update_progress(progress):
    import sys
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = '\r\n'
    #status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()
#
def compute_intv(v,t,p,p0=1e3):

    # some constants
    kappa = 2./7

    # pressure quantitites, p in hPa
    pp0 = (p0/p[np.newaxis,:,np.newaxis])**kappa
    dp  = np.gradient(p)[np.newaxis,:,np.newaxis]*100.
    # zonal means and eddy terms
    v_bar = v.mean(axis=-1)
    v = v - v_bar[:,:,:,np.newaxis] # v = v'
    t_bar = t.mean(axis=-1)*pp0 # t_bar = theta_bar
    t = v*(t*pp0[:,:,:,np.newaxis] - t_bar[:,:,:,np.newaxis]) # t = v'Th'

    dthdp = np.gradient(t_bar,1,dp,1)[1] # dthdp = d(theta_bar)/dp
    dthdp[dthdp==0] = np.NaN
    t = t/dthdp[:,:,:,np.newaxis] # t = v'Th'/(dTh_bar/dp)
    t_bar = t.mean(axis=-1) # t_bar = bar(v'Th'/(dTh_bar/dp))

    return v_bar,t_bar

##############################################################################################
def eof(X,n=1):
    """Principal Component Analysis / Empirical Orthogonal Functions / SVD
        
        Uses Singular Value Decomposition to find the dominant modes of variability.
        The field X can be reconstructed with Y = dot(EOF,PC) + X.mean(axis=time)
        
        INPUTS:
            X -- Field, shape (time x space).
            n -- Number of modes to extract
        OUTPUTS:
            EOF - Spatial modes of variability
            PC  - Temporal evolution of EOFs
            E   - Explained value of variability
            u   - spatial modes
            s   - variances
            v   - temporal modes
    """
    import numpy as np
    import scipy.signal as sg
    
    # find out which dimension is time
    #  assume that time is longer dimension
    #shpe = np.shape(X)
    #if shpe[0] > shpe[1]:
    #    X = X.T
    #    shpe = np.shape(X)
    # take out the time mean
    X = sg.detrend(X.T)
    # perform SVD - v is actually V.H in X = U*S*V.H
    u,s,v = np.linalg.svd(X, full_matrices=False)
    # now, u contains the spatial, and v the temporal structures
    # s contains the variances, with the same units as the input X
    # u.shape = (space, modes(space)), v.shape = (modes(space), time)

    # get the first n modes, in physical units
    #  we can either project the date onto  the principal component, F*V
    #  or multiply u*s. This is the same, as U*S*V.H*V = U*S
    EOF = np.dot(u[:,:n],diag(s)[:n,:n])
    # time evolution is in v
    PC  = v[:n,:]
    # EOF wants \lambda = the squares of the eigenvalues,
    #  but SVD yields \gamma = \sqrt{\lambda}
    s = s*s
    E   = s[:n]/sum(s)

    return EOF,PC,E,u[:,:n],sqrt(s[:n]),v.T[:,:n]


##############################################################################################
def ComputeVstar(inFileName, outFileName='same', temp='temp', vcomp='vcomp', lat='lat', pfull='pfull', time='time', p0=1e3, wkdir='./'):
    """Computes the residual meridional wind v* (as a function of time).
        
        INPUTS:
            inFileName  - filename of input file, relative to wkdir
            outFileName - filename of output file, 'none', or 'same'
            temp        - name of temperature field in inFile
            vcomp       - name of meridional velocity field in inFile
            lat         - name of latitude in inFile
            pfull       - name of pressure in inFile [hPa]
            time        - name of time field in inFile. Only needed if outFile used
            p0          - pressure basis to compute potential temperature [hPa]
            wkdir       - root directory, from with the paths inFile and outFile are taken
        OUTPUTS:
            vstar       - residual meridional wind, as a function of time
    """
    import netCDF4 as nc
    import numpy as np

    a0    = 6371000
    g     = 9.81

    # read input file
    print 'Reading data'
    update_progress(0)
    if outFileName == 'same':
        mode = 'a'
    else:
        mode = 'r'
    inFile = nc.Dataset(wkdir + inFileName, mode)
    t = inFile.variables[temp][:]
    update_progress(.25)
    v = inFile.variables[vcomp][:]
    update_progress(.50)
    l = inFile.variables[lat][:]
    update_progress(.75)
    p = inFile.variables[pfull][:]
    update_progress(1)


    v_bar,t_bar = compute_intv(v,t,p,p0) # t_bar = bar(v'Th'/(dTh_bar/dp))

    # Eulerian streamfunction
    dp  = np.gradient(p)[np.newaxis,:,np.newaxis]*100.
    vstar = v_bar - np.gradient(t_bar,1,dp,1)[1]

    return vstar
    


##############################################################################################
def GlobalAvg(lat,data,axis=-1,lim=20,mx=90):
    """Compute cosine weighted meridional average from lim to mx.

    INPUTS:
      lat  - latitude
      data - data to average N x latitude
      axis - axis designating latitude
      lim  - starting latitude to average
      mx   - stopping latitude to average
    OUTPUTS:
      integ- averaged data
    """
    from numpy import trapz,cos,prod,reshape,newaxis
    #get data into the correct shape
    tmp = AxRoll(data,axis)
    shpe= tmp.shape
    tmp = reshape(tmp,(shpe[0],prod(shpe[1:])))
    #cosine weighting
    J = find((lat>=lim)*(lat<=mx))
    coslat = cos(lat*pi/180.)
    coswgt = trapz(coslat[J],lat[J])
    tmp = trapz(tmp[J,:]*coslat[J,newaxis],lat[J],axis=0)/coswgt
    integ = reshape(tmp,shpe[1:])
    return integ

##############################################################################################
def GetWaves(x,y=[],wave=0,axis=-1,do_anomaly=False):
    from numpy import fft,squeeze
    x = AxRoll(x,axis)
    #compute anomalies
    if do_anomaly:
        x = GetAnomaly(x,0)
    if len(y) > 0:
        y = AxRoll(y,axis)
        if do_anomaly:
            y = GetAnomaly(y,0)
    #Fourier decompose
    x = fft.fft(x,axis=0)
    if wave == 0:
        xym = zeros_like(x)
    if len(y) > 0:
        y = fft.fft(y,axis=0)
        #Take out the waves - there's a magic factor of two
        nl  = x.shape[0]**2
        xym  = 2*real(x*y.conj())/nl
        if wave > 0:
            xym = xym[wave,:]/nl
    else:
        mask = zeros_like(x)
        if wave > 0:
            mask[wave,:] = 1
            xym = real(fft.ifft(x*mask,axis=0))
        else:
            xym = zeros_like(x)
            for m in range(x.shape[0]):
                mask[m,:] = 1
                xym[m,:] = real(fft.ifft(x*mask,axis=0)).mean(axis=0)
                mask[m,:] = 0
    xym = AxRoll(xym,axis,'i')
    return xym

##helper functions
def GetAnomaly(x,axis=-1):
    """Computes the anomaly of array x along dimension axis.

    INPUTS:
      x    - array to compute anomalies from
      axis - axis along dimension for anomalies
    OUTPUTS:
      x    - anomalous array
    """
    #bring axis to the front
    xt= AxRoll(x,axis)
    #compute anomalies
    xt = xt - xt.mean(axis=0)[newaxis,:]
    #bring axis back to where it was
    x = AxRoll(xt,0,'i')
    return x
#
def AxRoll(x,ax,start_mode=0):
    from numpy import rollaxis
    if isinstance(start_mode, basestring):
        mode = start_mode
    else:
        mode = 'f'
    #
    if ax < 0:
        n = len(x.shape) + ax
    else:
        n = ax
    #
    if mode is 'f':
        y = rollaxis(x,n,start_mode)
    elif mode is 'i':
        y = rollaxis(x,0,n+1)
    else:
        raise Exception("mode must be 'f' for forward or 'i' for inverse")
    return y
      
