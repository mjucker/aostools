#!/usr/bin/python
# Filename: climate.py
#
# Code by Martin Jucker, distributed under an GPLv3 License
# Any publication benefitting from this piece of code should cite
# Jucker, M 2014. Scientific Visualisation of Atmospheric Data with ParaView.
# Journal of Open Research Software 2(1):e4, DOI: http://dx.doi.org/10.5334/jors.al
#
# This file provides helper functions that can be useful as pre-viz-processing of files and data

############################################################################################
#
# compute climatologies

## helper function: check if string contained in upper level string
def CheckAny(set,string):
    for c in set:
        if c in string: return 1
    return 0

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
        outFile     name of the output file created
        writes outputfile with name depending on input file name and climatType
    """

    # need to read netCDF and of course do some math

    import netCDF4 as nc
    import numpy as np

    TimeStepRange = 1.01
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
    timeUnits = timeVar.units
    chck = CheckAny(('seconds','days','months'),timeUnits)
    if not chck:
        print 'Cannot understand units of time, which is: '+timeUnits
        timeUnits = raw_input('Please provide units [seconds,days,months] ')
    if 'seconds' in timeUnits:
        timeStepFact = 1./86400
        timeUnits = 'seconds'
    elif 'days' in timeUnits:
        timeStepFact = 1
        timeUnits = 'days'
    elif 'months' in timeUnits:
        timeStepFact = 30
        timeUnits = 'months'
    else:
        raise Exception("I don't understand the time unit: "+timeVar.units)
    if np.diff(timeVar).max()/np.diff(timeVar).min() > TimeStepRange:
        raise Exception("The time step is not constant, but varies between "+str(diff(timeVar).min())+" and "+str(diff(timeVar).max()))
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
        raise Exception("File "+wkdir+file+": Cannot deal with years that are not 360 or 365 days. This data has "+str(totalNumDays)+" days.")
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
        if len(varShape) == 0: continue
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
    return outFileName

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
            psi         - stream function, as a function of time
            psis        - residual stream function, as a function of time
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

    v_bar,t_bar = ComputeVertEddy(v,t,p,p0) # t_bar = bar(v'Th'/(dTh_bar/dp))

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


##############################################################################################
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
def ComputeVertEddy(v,t,p,p0=1e3,wave=-1):
    """ Computes the vertical eddy components of the residual circulation,
        bar(v'Theta'/Theta_p). Either in real space, or a given wave number.
        Dimensions must be time x pres x lat x lon
        
        INPUTS:
            v    - meridional wind
            t    - temperature
            p    - pressure coordinate
            p0   - reference pressure for potential temperature
            wave - wave number (if >=0)
        OUPUTS:
            v_bar - zonal mean meridional wind
            t_bar - zonal mean vertical eddy component <v'Theta'/Theta_p>
    """
    import numpy as np
    #
    # some constants
    kappa = 2./7
    #
    # pressure quantitites, p in hPa
    pp0 = (p0/p[np.newaxis,:,np.newaxis,np.newaxis])**kappa
    dp  = np.gradient(p)[np.newaxis,:,np.newaxis]*100.
    # convert to potential temperature
    t = t*pp0 # t = theta
    # zonal means
    v_bar = v.mean(axis=-1)
    t_bar = t.mean(axis=-1) # t_bar = theta_bar
    # prepare pressure derivative
    dthdp = np.gradient(t_bar,1,dp,1,edge_order=2)[1] # dthdp = d(theta_bar)/dp
    dthdp[dthdp==0] = np.NaN
    # now get wave component(s)
    if wave < 0:
        v = GetAnomaly(v)
        t = GetAnomaly(t)
        t = v*t
        t = t/dthdp[:,:,:,np.newaxis] # t = v'Th'/(dTh_bar/dp)
        t_bar = t.mean(axis=-1) # t_bar = bar(v'Th'/(dTh_bar/dp))
    else:
        t = GetWaves(v,t,wave=wave,do_anomaly=True) # t = v'Th'_{k=wave}
        t_bar = t/dthdp # t_bar = bar(v'Th'/(dTh_bar/dp))
    #
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
    #  we can either project the data onto the principal component, F*V
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
def ComputeAnnularMode(lat, pres, time, data, choice='z'):
    """Compute annular mode as in Geber et al, GRL 2008.
        This is basically the first PC, but normalized to unit variance and zero mean.
        
        INPUTS:
            lat    - latitude
            pres   - pressure
            time   - time
            data   - variable to compute EOF from. This is typically
                        geopotential or zonal wind.
                        Size time x pres x lat (ie zonal mean)
            choice - not essential, but used for sign convention.
                        If 'z', the sign is determined based on 70-80N.
                        Otherwise, 50-60N is used.
        OUTPUT:
            AM     - The annular mode, size time x pres
    """
    import numpy as np
    from matplotlib.mlab import find
    #
    AM = np.empty((len(time),len(pres)))
    AM[:] = np.nan
    j_tmp = find(lat > 20)
    coslat = np.cos(lat*np.pi/180.)
    negCos = (coslat < 0.)
    coslat[negCos] = 0.
    # weighting as in Gerber et al GRL 2008
    sqrtcoslat = np.sqrt(coslat[j_tmp])
    # try to get the sign right
    # first possibility
    if choice == 'z':
        minj = 70
        maxj = 80
        sig = -1
    else:
        minj = 50
        maxj = 60
        sig = 1
    jj = (lat[j_tmp] > minj)*(lat[j_tmp] < maxj)
    # second possibility
    #jj = abs(lat[j_tmp]-80).argmin()
    #sig = -1
    for k in range(len(pres)):
        # remove global mean
        globZ = GlobalAvg(lat,data[:,k,:],axis=-1)
        var = data[:,k,:] - globZ[:,np.newaxis]
        # area weighting: EOFs are ~variance, thus take sqrt(cos)
        var = var[:,j_tmp]*sqrtcoslat[np.newaxis,:]
        varNan = np.isnan(var)
        if np.sum(np.reshape(varNan,(np.size(varNan),)))==0:
            eof1,pc1,E,u,s,v = eof(var)
            # force the sign of PC
            pc1  = pc1*sig*np.sign(eof1[jj].mean())
            # force unit variance and zero mean
            AM[:,k] = (pc1-pc1.mean())/np.std(pc1)
    return AM

##############################################################################################
def ComputeVstar(data, temp='temp', vcomp='vcomp', pfull='pfull', wave=-1, p0=1e3, wkdir='./'):
    """Computes the residual meridional wind v* (as a function of time).
        
        INPUTS:
            data  - filename of input file, relative to wkdir, or dictionary with {T,v,pfull}
            temp  - name of temperature field in data
            vcomp - name of meridional velocity field in data
            pfull - name of pressure in inFile [hPa]
            wave  - decompose into given wave number contribution if wave>=0
            p0    - pressure basis to compute potential temperature [hPa]
            wkdir - root directory, from with the path to the input file is taken
        OUTPUTS:
            vstar       - residual meridional wind, as a function of time
    """
    import netCDF4 as nc
    import numpy as np

    a0    = 6371000
    g     = 9.81

    # read input file
    if isinstance(data,str):
        print 'Reading data'
        update_progress(0)
        #
        inFile = nc.Dataset(wkdir + data, 'r')
        t = inFile.variables[temp][:]
        update_progress(.45)
        v = inFile.variables[vcomp][:]
        update_progress(.90)
        p = inFile.variables[pfull][:]
        update_progress(1)
        inFile.close()
        #
        v_bar,t_bar = ComputeVertEddy(v,t,p,p0,wave=wave)
    else:
        p = data[pfull]
        v_bar,t_bar = ComputeVertEddy(data[vcomp],data[temp],p,p0,wave=wave)
    # t_bar = bar(v'Th'/(dTh_bar/dp))
    #
    dp  = np.gradient(p)[np.newaxis,:,np.newaxis]*100.
    vstar = v_bar - np.gradient(t_bar,1,dp,1)[1]

    return vstar
    

##############################################################################################
def ComputeWstar(data, slice='all', omega='omega', temp='temp', vcomp='vcomp', pfull='pfull', lat='lat', wave=[-1], p0=1e3):
    """Computes the residual upwelling w* as a function of time.
    	input dimensions must be time x pres x lat x lon.
    	output is either space-time (wave<0, dimensions time x pres x lat)
         or space-time-wave (dimensions wave x time x pres x lat).
        
        INPUTS:
            data  - filename of input file, relative to wkdir, or dictionary with (w,T,v,pfull)
            slice - time slice to work with (large memory requirements). Array [start,stop] or 'all'
            omega - name of pressure velocity field in data
            temp  - name of temperature field in data
            vcomp - name of meridional velocity field in data
            pfull - name of pressure in data [hPa]
            lat   - name of latitude in data [deg]
            wave  - decompose into given wave number contribution(s) if 
            		 len(wave)=1 and wave>=0, or len(wave)>1
            p0    - pressure basis to compute potential temperature [hPa]
        OUTPUTS:
            w_bar - zonal mean Eulerian pressure velocity
            wstar - residual pressure velocity, as a function of time [and waves]
    """
    import netCDF4 as nc
    import numpy as np

    a0    = 6371000.
    
    # read input file
    if isinstance(data,str):
        print 'Reading data'
        #
        #if outFile == 'same': outFile = data
        inFile = nc.Dataset(data, 'r')
        if slice == 'all':
        	slice=[0,inFile.variables[omega][:].shape[0]]
        data = {}
        data[omega] = inFile.variables[omega][slice[0]:slice[1],:]
        data[temp] = inFile.variables[temp][slice[0]:slice[1],:]
        data[vcomp] = inFile.variables[vcomp][slice[0]:slice[1],:]
        data[pfull] = inFile.variables[pfull][:]
        data[lat] = inFile.variables[lat][:]
        inFile.close()
    # spherical geometry
    pilat = data[lat]*np.pi/180.
    coslat = np.cos(pilat)[np.newaxis,np.newaxis,:]
    R = a0*coslat[np.newaxis,:]
    R = 1./R
    dphi = np.gradient(pilat)[np.newaxis,np.newaxis,:]
    # compute thickness weighted meridional heat flux
    shpe = data[omega].shape[:-1]
    vt_bar = zeros((len(wave),)+shpe)
    for w in range(len(wave)):
    	#  w_bar is actually v_bar, but we don't need that
    	w_bar,vt_bar[w,:] = ComputeVertEddy(data[vcomp],data[temp],data[pfull],p0,wave=wave[w])
    	# weigh v'T' by cos\phi
    	vt_bar[w,:] = vt_bar[w,:]*coslat
    	# get the meridional derivative
    	vt_bar[w,:] = np.gradient(vt_bar[w,:],1,1,dphi,edge_order=2)[-1]
    # compute zonal mean upwelling
    w_bar = data[omega].mean(axis=-1)
    # put it all together
    if len(wave)==1:
    	return w_bar,np.squeeze(R*vt_bar)
    else:
    	return w_bar,R*vt_bar
##############################################################################################
def ComputeEPfluxDiv(lat,pres,upvp,vptp,tbar,ubar=None,upwp=None):
    """ Compute the EP-flux vectors and divergence terms.

        The vectors are normalized to be plotted in cartesian (linear)
        coordinates, i.e. do not include the geometric factor a*cos\phi.
        Thus, ep1 is in [m2/s2], and ep2 in [hPa*m/s2].
        The divergence is in units of m/s/day, and therefore represents
        the deceleration of the zonal wind. This is actually the quantity
        1/(acos\phi)*div(F).

    INPUTS:
      lat    - latitude [degrees]
      pres   - pressure [hPa]
      upvp   - [u'v'], shape(time,p,lat)
      vptp   - [v'T'], shape(time,p,lat)
      tbar   - [T]   , shape(time,p,lat)
      ubar   - [u]   , optional, shape(time,p,lat)
      upwp   - [u'w'], optional, shape(time,p,lat)
    OUTPUTS:
      ep1  - meridional EP-flux component, scaled to plot in cartesian [m2/s2]
      ep2  - vertical   EP-flux component, scaled to plot in cartesian [hPa*m/s2]
      div1 - horizontal EP-flux divergence, divided by acos\phi [m/s/d]
      div2 - horizontal EP-flux divergence , divided by acos\phi [m/s/d]
    """
    from numpy import cos,tan,mat,newaxis,reshape,mat,gradient
    # some constants
    R     = 287.04
    cp    = 1004
    kappa = R/cp
    p0    = 1000
    Omega = 2*pi/(24*3600)
    a0    = 6.378e6
    # geometry
    pilat = lat*pi/180
    dphi  = gradient(pilat)[newaxis,newaxis,:]
    coslat= cos(pilat)[newaxis,newaxis,:]
    sinlat= sin(pilat)[newaxis,newaxis,:]
    tanlat= tan(pilat)[newaxis,newaxis,:]
    f     = 2*Omega*sinlat
    pp0  = (p0/pres[newaxis,:,newaxis])**kappa
    dp    = gradient(pres)[newaxis,:,newaxis]
    #
    # absolute vorticity
    if ubar is None:
        fhat = 0.
    else:
        fhat = gradient(ubar*coslat,1,1,dphi,edge_order=2)[-1]/coslat/a0
    fhat = f - fhat
    #
    ## compute thickness weighted heat flux
    # dtdp = [theta]_p:
    theta = tbar*pp0
    dtdp  = gradient(theta,1,dp,1,edge_order=2)[1]
    # vptp = [v'theta']
    vptp = vptp*pp0  
    #
    ## compute the horizontal component
    if ubar is None:
        shear = 0.
    else:
        shear = gradient(ubar,1,dp,1,edge_order=2)[1]
    ep1_cart = -upvp + shear*vptp/dtdp
    #
    ## compute vertical component of EP flux.
    ## at first, keep it in Cartesian coordinates, ie ep2_cart = f [v'theta']_bar / [theta]_p + ...
    #
    ep2_cart = fhat*vptp/dtdp
    if upwp is not None:
        ep2_cart -= upwp
    #
    #
    #
    # div1 = 1/a0 d/d phi [-u'v'] - 2/a0 tan(phi)[-u'v']
    # div2 = d/dp (f[v'theta']/[theta]_p)
    div1 = gradient(ep1_cart,1,1,dphi,edge_order=2)[-1] - 2*tanlat*ep1_cart
    div1 = div1/a0
    #
    div2 = gradient(ep2_cart,1,dp,1,edge_order=2)[1]
    #
    # convert to m/s/day
    div1 = div1*86400
    div2 = div2*86400
    #
    return ep1_cart,ep2_cart,div1,div2

##############################################################################################
def GlobalAvg(lat,data,axis=-1,lim=20,mx=90,cosp=1):
    """Compute cosine weighted meridional average from lim to mx.

    INPUTS:
      lat  - latitude
      data - data to average N x latitude
      axis - axis designating latitude
      lim  - starting latitude to average
      mx   - stopping latitude to average
      cosp - power of cosine weighting
    OUTPUTS:
      integ- averaged data
    """
    from numpy import trapz,cos,prod,reshape,newaxis,pi
    #get data into the correct shape
    tmp = AxRoll(data,axis)
    shpe= tmp.shape
    tmp = reshape(tmp,(shpe[0],prod(shpe[1:])))
    #cosine weighting
    J = find((lat>=lim)*(lat<=mx))
    coslat = cos(lat*pi/180.)**cosp
    coswgt = trapz(coslat[J],lat[J])
    tmp = trapz(tmp[J,:]*coslat[J][:,newaxis],lat[J],axis=0)/coswgt
    integ = reshape(tmp,shpe[1:])
    return integ

##############################################################################################
def GetWaves(x,y=[],wave=-1,axis=-1,do_anomaly=False):
	"""Get Fourier mode decomposition of x, or x*y
	
	INPUTS:
		x          - the array to decompose
		y          - second array if wanted
		wave       - which mode to extract. all if -1
		axis       - along which axis of x (and y) to decompose
		do_anomaly - decompose from anomalies or full data
	OUTPUTS:
		xym        - data in Fourier space
	"""
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
	if wave < 0:
		xym = zeros_like(x)
	if len(y) > 0:
		y = fft.fft(y,axis=0)
        #Take out the waves - there's a magic factor of two
		nl  = x.shape[0]**2
		xym  = 2*real(x*y.conj())/nl
		if wave >= 0:
			xym = xym[wave,:][newaxis,:]
	else:
		mask = zeros_like(x)
		if wave >= 0:
			mask[wave,:] = 1
			xym = real(fft.ifft(x*mask,axis=0))
		else:
			for m in range(x.shape[0]):
				mask[m,:] = 1
				xym[m,:] = real(fft.ifft(x*mask,axis=0)).mean(axis=0)
				mask[m,:] = 0
	xym = AxRoll(xym,axis,'i')
	return squeeze(xym)

##helper functions
def GetAnomaly(x,axis=-1):
    """Computes the anomaly of array x along dimension axis.

    INPUTS:
      x    - array to compute anomalies from
      axis - axis along dimension for anomalies
    OUTPUTS:
      x    - anomalous array
    """
    from numpy import newaxis
    #bring axis to the front
    xt= AxRoll(x,axis)
    #compute anomalies
    xt = xt - xt.mean(axis=0)[newaxis,:]
    #bring axis back to where it was
    x = AxRoll(xt,axis,'i')
    return x
#
def AxRoll(x,ax,start_mode=0):
    """Re-arrange array x so that axis 'ax' is first dimension.
        Undo this if start_mode=='i'
    """
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


#######################################################
def Meters2Coord(data,mode='m2lat',coord=[],axis=-1):
    """Convert value (probably vector component) of one unit
        into another one.
        
        INPUTS:
        data  - data to convert
        mode  - 'm2lat', 'm2lon', 'm2hPa', and all inverses
        coord - values of latitude [degrees] or pressure [hPa]
        axis  - axis of data which needs modification
        OUTPUTS:
        out   - converted from data
        """
    from numpy import cos,pi
    # constants
    a0    = 6.378e6
    ps    = 1e3
    H     = 7e3
    # geometric quantitites
    rad2deg = 180/pi
    coslat  = cos(lat/rad2deg)
    cosm1   = 1/coslat
    gemfac  = rad2deg/a0
    #
    if mode is 'm2lat':
        out = data*gemfac
    elif mode is 'lat2m':
        out = data/gemfac
    elif mode in ['m2lon','lon2m','m2hPa','hPa2m']:
        tmp = AxRoll(data,axis)
        out = zeros_like(tmp)
    else:
        raise ValueError("mode not recognized")
    if mode is 'm2lon':
        for l in range(out.shape[0]):
            out[l,:] = tmp[l,:]*cosm1
        out = out*gemfac
        out = AxRoll(out,axis,'i')
    elif mode is 'lon2m':
        for l in range(out.shape[0]):
            out[l,:] = tmp[l,:]*coslat
        out = out/gemfac
    elif mode is 'm2hPa':
        for p in range(out.shape[0]):
            out[p,:] = -coord[p]*tmp[p,:]
        out = out/H
        out = AxRoll(out,axis,'i')
    elif mode is 'hPa2m':
        for p in range(out.shape[0]):
            out[p,:] = -coord[p]/tmp[p,:]
        out[tmp==0] = NaN
        out = out*H
        out = AxRoll(out,axis,'i')
    #
    return out

##############################################################################################
def ComputeBaroclinicity(lat, tempIn, hemi='both', minLat=20, maxLat=60, pres=None, minPres=250):
    """Compute the meridional temperature gradient, integrated between minLat and maxLat.
    	Thus, baroclinicity is defined as the difference in temperature between minLat and maxLat,
    	optionally integrated from the surface to minPres.
    	
        INPUTS:
        lat     - latitude [degrees] 
        tempIn  - temperature, shape time x pressure x lat
        hemi    - hemisphere to consider. 'both','N','S'
        minLat  - latitude closest to equator
        maxLat  - latitdue closest to pole
        pres    - pressure [hPa]
        minPres - top of pressure averaging from surface.
        			only used if pres is not None.
        OUTPUT:
        dT      - dictionary with options ['S'] and ['N'] if hemi='both'
        			otherwise array of length len(time)
    """
    import numpy as np

    numTimeSteps = tempIn.shape[0]
    # get important meridional grid points
    N = find( (lat>= minLat)*(lat<= maxLat) )
    #  make sure first index is closer to pole
    if lat[N[0]] < lat[N[-1]]:
        latN = [N[0],N[-1]]
    else:
        latN = [N[-1],N[0]]
    S = find( (lat<=-minLat)*(lat>=-maxLat) )
    #  make sure first index is closer to pole
    if lat[S[0]] > lat[S[-1]]:
        latS = [S[0],S[-1]]
    else:
        latS = [S[-1],S[0]]
    # pressure levels
    if pres is None:
        temp = tempIn[:,np.newaxis,:]
        K = 0
    else:
        temp = tempIn
        K = find(pres >= minPres)
    #
    T = dict()
    for H in ['S','N']:
        T[H] = np.zeros((temp.shape[0],len(K),2))
    for l in range(2):
        T['S'][:,:,l] = temp[:,K,latS[l]]
        T['N'][:,:,l] = temp[:,K,latN[l]]
    #
    dT = dict()
    for H in ['S','N']:
        if pres is None:
            tmp = T[H]
        else:
            tmp = trapz(T[H],pres[K],axis=1)/trapz(pres[K])
        # zero index is closer to tropics
        #  define bariclinicity as -dyT
        dT[H] = tmp[:,0] - tmp[:,1]
    if hemi == 'both':
        return dT
    else:
        return dT[hemi]


#######################################################
def SymmetricColorbar(fig, obj, zero=0):
	"""Make colorbar symmetric with respect to zero.
		Note: this does not update the colobar limits, but adjusts
		the colormap such that the node is around zero.
		
		INPUTS:
		fig  - figure object to attach colobar to
		obj  - color plot [e.g. obj=contourf()]
		zero - node to assign color of 0
	"""
	cb = fig.colorbar(obj)
	(cmn,cmx) = cb.get_clim()
	cmnp = cmn-zero
	cmxp = cmx-zero
	c0 = min(cmnp,-cmxp) + zero
	c1 = max(cmxp,-cmnp) + zero
	obj.set_clim(c0,c1)






