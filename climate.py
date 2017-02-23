#!/usr/bin/python
# Filename: climate.py
#
# Code by Martin Jucker, distributed under an GPLv3 License
#
# This file provides helper functions that can be useful as pre-viz-processing of files and data

############################################################################################
#
# compute climatologies

## helper function: check if string contained in list (set) of strings
def CheckAny(string,set):
    for c in set:
        if c in string: return True
    return False

## helper function: return the day of the year instead of full date
def FindDayOfYear(dateStruc,dateUnits,calendar):
    import netcdftime as nct
    nDays = len(dateStruc)
    t = nct.utime(dateUnits,calendar=calendar)
    dateLoc = np.zeros_like(dateStruc)
    for d in range(nDays):
        dateLoc[d] = nct.datetime(1,dateStruc[d].month,dateStruc[d].day)
    dayOfYear = t.date2num(dateLoc)
    return dayOfYear

## compute climatologies
def ComputeClimate(file, climatType, wkdir='/', timeDim='time',cal=None):
    """Compute climatologies from netCDF files.

        ComputeClimate(file,climatType,wkdir='/',timeDim='time')

        Inputs:
        file        file name, relative path from wkdir
        climatType  'daily', 'monthly', 'annual', 'DJF', 'JJA', or any
                    combination of months according to two-letter code
                    Ja Fe Ma Ap My Jn Jl Au Se Oc No De
        wkdir       working directory, in which 'file' must be, and to which the output
                    is written
        timeDim     name of the time dimension in the netcdf file
        cal         calendar, if other than within the netcdf file
        Outputs:
        outFile     name of the output file created
        writes outputfile with name depending on input file name and climatType
    """

    # need to read netCDF and of course do some math
    import netCDF4 as nc
    import numpy as np
    import os

    if climatType == 'DJF':
        climType = 'DeJaFe'
    elif climatType == 'JJA':
        climType = 'JuJlAu'
    elif climatType == 'annual':
        climType = 'JaFeMaApMyJnJlAuSeOcNoDe'
    else:
        climType = climatType
    monthList=['Ja','Fe','Ma','Ap','My','Jn','Jl','Au','Se','Oc','No','De']
    calendar_types = ['standard', 'gregorian', 'proleptic_gregorian', 'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day']

    if wkdir[-1] != '/': wkdir += '/'
    if os.path.isfile(wkdir+file):
        ncFile = nc.Dataset(wkdir+file,'r+')
    else:
        raise IOError(wkdir+file+' does not exist')

    time    = ncFile.variables[timeDim][:]
    numTimeSteps = len(time)
    timeVar = ncFile.variables[timeDim]
    # check the time units
    timeUnits = timeVar.units
    chck = CheckAny(timeUnits,('seconds','days','months'))
    if not chck:
        print 'Cannot understand units of time, which is: '+timeUnits
        newUnits = raw_input('Please provide units [seconds,days,months] ')
        if newUnits not in ["seconds","days","months"]:
            raise ValueError('units must be seconds, days, or months')
        unitSplit = timeUnits.split()
        unitSplit[0] = newUnits
        timeUnits = ' '.join(unitSplit)
    timeStep = np.diff(timeVar).mean()
    print 'The time dimension is in units of',timeUnits,', with a mean time step of',timeStep,'days'
    # check the calendar type
    getCal = False
    if cal:
        timeCal = cal
    else:
        try:
            timeCal = str(timeVar.calendar)
            if not CheckAny(timeCal,calendar_types):
                print 'Cannot understand the calendar type, which is: '+timeCal
                timeCal = raw_input('Please provide a calendar type from the list '+str(calendar_types)+' ')
                timeVar.calendar = timeCal
        except:
            timeCal = raw_input('Please provide a calendar type from the list '+str(calendar_types)+' ')
    if timeCal not in calendar_types:
        raise ValueError('calender must be in '+str(calendar_types))
    else:
        print 'Calendar type '+timeCal
    #
    # split everything into years,months,days
    date = nc.num2date(time,timeUnits,timeCal)
    days = np.zeros(len(date),)
    monthsI = np.zeros_like(days)
    monthsS = []
    years  = np.zeros_like(days)
    for d in range(len(date)):
        days[d] = date[d].day
        monthsI[d] = date[d].month
        monthsS.append(monthList[date[d].month-1])
        years[d] = date[d].year
    # Now, need to know about the type of climatology we want.
    #
    if climType == 'daily':
        dayOfYear = FindDayOfYear(date,timeUnits,timeCal)
        climTimeDim = np.sort(np.unique(dayOfYear))
        climTimeVar = dayOfYear
    elif climType == 'monthly':
        climTimeDim = np.sort(np.unique(monthsI)) - 1
        climTimeVar = monthsI - 1
    else:
        climTimeVar = np.zeros_like(days)
        for m in range(len(climType)/2):
            thisMonth = climType[m*2:m*2+2]
            indices = [i for i, x in enumerate(monthsS) if x == thisMonth]
            climTimeVar[indices] = 1

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
                if not 'edges' in att:
                    outVar.setncattr(att,inVar.getncattr(att))
        elif climType == 'daily' or climType == 'monthly':
            nTime = len(climTimeDim)
            if climType == 'daily':
                units = 'days'
            else:
                units = 'months'
            dTime = climTimeDim
            outDim = outFile.createDimension(dim,nTime)
            timeValue = dTime
            outVar = outFile.createVariable(dim,str(ncFile.variables[dim].dtype),(dim,))
            outVar[:] = timeValue
            outVar.setncattr('long_name','climatological ' + units[:-1] + ' of year')
            outVar.setncattr('units',units + ' since 0001-01-01 00:00:00')
            outVar.setncattr('calendar',timeCal)
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
                tmpAvg = tmpVar[climTimeVar>0,:].mean(axis=0)
            else:
                outVar = outFile.createVariable(var,str(ncFile.variables[var].dtype),ncFile.variables[var].dimensions    )
                avgShape = []
                avgShape.append(nTime)
                for t in range(len(np.shape(outVar))-1):
                    avgShape.append(np.shape(outVar)[t+1])
                tmpAvg = np.zeros(avgShape)
                for t in range(nTime):
                    includeSteps = climTimeVar == climTimeDim[t]
                    tmpAvg[t,:] = tmpVar[includeSteps,:].mean(axis=0)
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
# get the saturation mixing ration according to Clausius-Clapeyron

# helper function: re-arrange array dimensions
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


def ComputeSaturationMixingRatio(T, p, pDim):
    """Computes the saturation water vapor mixing ratio according to Clausius-Clapeyron

        INPUTS:
            T    - temperature in Kelvin, any size
            p    - pressure in hPa/mbar, must be one dimension of T
            pDim - index of dimension corresponding to p
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

    # make sure we are operating along the pressure axis
    T = AxRoll(T,pDim)

    # pressure is assumed in hPa: convert to Pa
    p = p*100
    # compute saturation pressure
    esat = ES0*np.exp(HLV*(1./Tf - 1./T)/Rv)
    qsat = np.zeros_like(esat)
    # finally, compute saturation mixing ratio from pressure
    for k in range(len(p)):
        qsat[k,:] = Rd/Rv*esat[k,:]/(p[k]-esat[k,:])
    return AxRoll(qsat,pDim,'i')


##############################################################################################
def ComputeRelativeHumidity(inFile, outFile='none', temp='temp', sphum='sphum', pfull='pfull'):
    """Computes relative humidity from temperature and specific humidity.

        File inFile is assumed to contain both temperature and specific humidity.
        Relative humidity is either output of the function, or written to the file outFile.

        Inputs:
            inFile    Name of the file (full path)
                        containing temperature and moisture
            outFile   Name of the output file containing specific humidity.
                        No output file is created if outFile='none'
            temp      Name of the temperature variable inside inFile
            sphum     Name of specific humidity variable inside inFile
            pfull     Name of full level pressure [hPa] inside inFile
    """

    import netCDF4 as nc
    import numpy as np
    # relative humidity is then q/qsat*100[->%]

    # read input file
    inFile = nc.Dataset(inFile, 'r')
    t = inFile.variables[temp][:]
    q = inFile.variables[sphum][:]
    p = inFile.variables[pfull][:]

    # compute saturation mixing ratio
    qsat = ComputeSaturationMixingRatio(t, p)

    #write output file
    if outFile is not 'none':
        outFile = nc.Dataset(inFile[0:-3]+'_out.nc','w')
        for dim in ncFile.dimensions:
            outDim = outFile.createDimension(dim,len(ncFile.dimensions[dim]))
            inVar = ncFile.variables[dim]
            outVar = outFile.createVariable(dim, str(ncFile.variables[dim].dtype),(dim,))
            for att in inVar.ncattrs():
                outVar.setncattr(att,inVar.getncattr(att))
            outVar[:] = inVar[:]
        outVar = outFile.createVariable('rh', 'f4', ncFile.variables[temp].dimensions)
        outVar[:] = q/qsat*1.e2
    return q/qsat*1.e2


##############################################################################################
def ComputePsi(data, outFileName='none', temp='temp', vcomp='vcomp', lat='lat', pfull='pfull', time='time', p0=1e3):
    """Computes the residual stream function \Psi* (as a function of time).

        INPUTS:
            data        - filename of input file or dictionary with temp,vcomp,lat,pfull
            outFileName - filename of output file, 'none', or 'same'
            temp        - name of temperature field in inFile
            vcomp       - name of meridional velocity field in inFile
            lat         - name of latitude in inFile
            pfull       - name of pressure in inFile [hPa]
            time        - name of time field in inFile. Only needed if outFile used
            p0          - pressure basis to compute potential temperature [hPa]
        OUTPUTS:
            psi         - stream function, as a function of time
            psis        - residual stream function, as a function of time
    """
    import netCDF4 as nc
    import numpy as np
    from scipy.integrate import cumtrapz
    import os

    # some constants
    kappa = 2./7
    a0    = 6371000
    g     = 9.81

    if isinstance(data,str):
        # check if file exists
        if not os.path.isfile(data):
            raise IOError('File '+data+' does not exist')
        # read input file
        print 'Reading data'
        update_progress(0)
        if outFileName == 'same':
            mode = 'a'
        else:
            mode = 'r'
        inFile = nc.Dataset(data, mode)
        t = inFile.variables[temp][:]
        update_progress(.45)
        v = inFile.variables[vcomp][:]
        update_progress(.90)
        l = inFile.variables[lat][:]
        update_progress(.95)
        p = inFile.variables[pfull][:]
        update_progress(1)
    else:
        t = data[temp]
        v = data[vcomp]
        l = data[lat]
        p = data[pfull]
        data = []
    p = p*100 # [Pa]
    #
    ## compute psi

    v_bar,t_bar = ComputeVertEddy(v,t,p,p0) # t_bar = bar(v'Th'/(dTh_bar/dp))

    # Eulerian streamfunction
    psi = cumtrapz(v_bar,x=p,axis=1,initial=0) # [m.Pa/s]
    v_bar=v=t=[]


    ## compute psi* = psi - bar(v'Th'/(dTh_bar/dp))
    psis = psi - t_bar
    t_bar = []
    psi = 2*np.pi*a0/g*psi *np.cos(l[np.newaxis,np.newaxis,:]*np.pi/180.) #[kg/s]
    psis= 2*np.pi*a0/g*psis*np.cos(l[np.newaxis,np.newaxis,:]*np.pi/180.) #[kg/s]

    ## write outputfile
    if outFileName is not 'none':
        print 'Writing file  '+outFileName
        if outFileName is not 'same':
            outFile = nc.Dataset(outFileName,'w')
            for dim in inFile.dimensions:
                if dim in [time,pfull,lat]:
                    outDim = outFile.createDimension(dim,len(inFile.dimensions[dim]))
                    inVar = inFile.variables[dim]
                    outVar = outFile.createVariable(dim, str(inFile.variables[dim].dtype),(dim,))
                    for att in inVar.ncattrs():
                        if att != '_FillValue': #no fill value in dimensions!
                            outVar.setncattr(att,inVar.getncattr(att))
                    outVar[:] = inVar[:]
        else:
            outFile = inFile
        outVar = outFile.createVariable('psi', 'f4', (time,pfull,lat,))
        outVar[:] = psi
        outVar = outFile.createVariable('psi_star', 'f4', (time,pfull,lat,))
        outVar[:] = psis
        outFile.close()
        print 'Done writing file '+outFileName
        if outFileName is not 'same':
            inFile.close()
    return psi,psis


##############################################################################################
## helper functions
def update_progress(progress,barLength=10):
    import sys
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
    text = "\r[{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), int(progress*100), status)
    sys.stdout.write(text)
    sys.stdout.flush()
#
def ComputeVertEddy(v,t,p,p0=1e3,wave=-1):
    """ Computes the vertical eddy components of the residual circulation,
        bar(v'Theta'/Theta_p). Either in real space, or a given wave number.
        Dimensions must be time x pres x lat x lon.
        Output dimensions are: time x pres x lat
        Output units are [v_bar] = [v], [t_bar] = [v*p]

        INPUTS:
            v    - meridional wind
            t    - temperature
            p    - pressure coordinate
            p0   - reference pressure for potential temperature
            wave - wave number (if >=0)
        OUPUTS:
            v_bar - zonal mean meridional wind [v]
            t_bar - zonal mean vertical eddy component <v'Theta'/Theta_p> [v*p]
    """
    import numpy as np
    #
    # some constants
    kappa = 2./7
    #
    # pressure quantitites
    pp0 = (p0/p[np.newaxis,:,np.newaxis,np.newaxis])**kappa
    dp  = np.gradient(p)[np.newaxis,:,np.newaxis]
    # convert to potential temperature
    t = t*pp0 # t = theta
    # zonal means
    v_bar = np.nanmean(v,axis=-1)
    t_bar = np.nanmean(t,axis=-1) # t_bar = theta_bar
    # prepare pressure derivative
    dthdp = np.gradient(t_bar,1,dp,1,edge_order=2)[1] # dthdp = d(theta_bar)/dp
    dthdp[dthdp==0] = np.NaN
    # time mean of d(theta_bar)/dp
    dthdp = np.nanmean(dthdp,axis=0)[np.newaxis,:]
    # now get wave component
    if wave < 0:
        v = GetAnomaly(v) # v = v'
        t = GetAnomaly(t) # t = t'
        t = np.nanmean(v*t,axis=-1) # t = bar(v'Th')
        t_bar = t/dthdp # t_bar = bar(v'Th')/(dTh_bar/dp)
    else:
        t = GetWaves(v,t,wave=wave,do_anomaly=True) # t = bar(v'Th'_{k=wave})
        t_bar = t/dthdp # t_bar = bar(v'Th')/(dTh_bar/dp)
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
    #
    AM = np.empty((len(time),len(pres)))
    AM[:] = np.nan
    j_tmp = np.where(lat > 20)[0]
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
def ComputeVstar(data, temp='temp', vcomp='vcomp', pfull='pfull', wave=-1, p0=1e3):
    """Computes the residual meridional wind v* (as a function of time).

        INPUTS:
            data  - filename of input file, relative to wkdir, or dictionary with {T,v,pfull}
            temp  - name of temperature field in data
            vcomp - name of meridional velocity field in data
            pfull - name of pressure in inFile [hPa]
            wave  - decompose into given wave number contribution if wave>=0
            p0    - pressure basis to compute potential temperature [hPa]
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
        inFile = nc.Dataset(data, 'r')
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
    dp  = np.gradient(p)[np.newaxis,:,np.newaxis]
    vstar = v_bar - np.gradient(t_bar,1,dp,1,edge_order=2)[1]

    return vstar


##############################################################################################
def ComputeWstar(data, slice='all', omega='omega', temp='temp', vcomp='vcomp', pfull='pfull', lat='lat', wave=[-1], p0=1e3):
    """Computes the residual upwelling w* as a function of time.

    	Input dimensions must be time x pres x lat x lon.
    	Output is either space-time (wave<0, dimensions time x pres x lat)
         or space-time-wave (dimensions wave x time x pres x lat).
        Output units are hPa/s, and the units of omega are expected to be hPa/s.

        INPUTS:
            data  - filename of input file, or dictionary with (w,T,v,pfull,lat)
            slice - time slice to work with (large memory requirements). Array [start,stop] or 'all'
            omega - name of pressure velocity field in data [hPa/s]
            temp  - name of temperature field in data
            vcomp - name of meridional velocity field in data
            pfull - name of pressure in data [hPa]
            lat   - name of latitude in data [deg]
            wave  - decompose into given wave number contribution(s) if
            		 len(wave)=1 and wave>=0, or len(wave)>1
            p0    - pressure basis to compute potential temperature [hPa]
        OUTPUTS:
            residual pressure velocity, time x pfull x lat [and waves] [hPa/s]
    """
    import netCDF4 as nc
    import numpy as np

    a0    = 6371000.

    # read input file
    if isinstance(data,str):
        inFile = nc.Dataset(data, 'r')
        if slice == 'all':
        	slice=[0,inFile.variables[omega][:].shape[0]]
        data = {}
        data[omega] = inFile.variables[omega][slice[0]:slice[1],:]*0.01 # [hPa/s]
        data[temp] = inFile.variables[temp][slice[0]:slice[1],:]
        data[vcomp] = inFile.variables[vcomp][slice[0]:slice[1],:]
        data[pfull] = inFile.variables[pfull][:] # [hPa]
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
    w_bar = np.nanmean(data[omega],axis=-1)
    # put it all together
    if len(wave)==1:
    	return w_bar + np.squeeze(R*vt_bar)
    else:
    	return w_bar + R*vt_bar
##############################################################################################
def ComputeEPfluxDiv(lat,pres,u,v,t,w=None,do_ubar=False,wave=-1):
    """ Compute the EP-flux vectors and divergence terms.

        The vectors are normalized to be plotted in cartesian (linear)
        coordinates, i.e. do not include the geometric factor a*cos\phi.
        Thus, ep1 is in [m2/s2], and ep2 in [hPa*m/s2].
        The divergence is in units of m/s/day, and therefore represents
        the deceleration of the zonal wind. This is actually the quantity
        1/(acos\phi)*div(F).

    INPUTS:
      lat  - latitude [degrees]
      pres - pressure [hPa]
      u    - zonal wind, shape(time,p,lat,lon) [m/s]
      v    - meridional wind, shape(time,p,lat,lon) [m/s]
      t    - temperature, shape(time,p,lat,lon) [K]
      w    - pressure velocity, optional, shape(time,p,lat,lon) [hPa/s]
      do_ubar - compute shear and vorticity correction? optional
      wave - only include this wave number. all if <0. optional
    OUTPUTS:
      ep1  - meridional EP-flux component, scaled to plot in cartesian [m2/s2]
      ep2  - vertical   EP-flux component, scaled to plot in cartesian [hPa*m/s2]
      div1 - horizontal EP-flux divergence, divided by acos\phi [m/s/d]
      div2 - horizontal EP-flux divergence , divided by acos\phi [m/s/d]
    """
    from numpy import pi,cos,sin,newaxis,gradient
    # some constants
    Rd    = 287.04
    cp    = 1004
    kappa = Rd/cp
    p0    = 1000
    Omega = 2*pi/(24*3600.) # [1/s]
    a0    = 6.371e6
    # geometry
    pilat = lat*pi/180
    dphi  = gradient(pilat)[newaxis,newaxis,:]
    coslat= cos(pilat)[newaxis,newaxis,:]
    sinlat= sin(pilat)[newaxis,newaxis,:]
    R     = 1./(a0*coslat)
    f     = 2*Omega*sinlat
    pp0  = (p0/pres[newaxis,:,newaxis])**kappa
    dp    = gradient(pres)[newaxis,:,newaxis]
    #
    # absolute vorticity
    if do_ubar:
        ubar = np.nanmean(u,axis=-1)
        fhat = R*gradient(ubar*coslat,1,1,dphi,edge_order=2)[-1]
    else:
        fhat = 0.
    fhat = f - fhat # [1/s]
    #
    ## compute thickness weighted heat flux [m.hPa/s]
    vbar,vertEddy = ComputeVertEddy(v,t,pres,p0,wave) # vertEddy = bar(v'Th'/(dTh_bar/dp))
    #
    ## get zonal anomalies
    u = GetAnomaly(u)
    v = GetAnomaly(v)
    if wave<0:
        upvp = np.nanmean(u*v,axis=-1)
    else:
        upvp = GetWaves(u,v,wave=wave)
    #
    ## compute the horizontal component
    if do_ubar:
        shear = gradient(ubar,1,dp,1,edge_order=2)[1] # [m/s.hPa]
    else:
        shear = 0.
    ep1_cart = -upvp + shear*vertEddy # [m2/s2 + m/s.hPa*m.hPa/s] = [m2/s2]
    #
    ## compute vertical component of EP flux.
    ## at first, keep it in Cartesian coordinates, ie ep2_cart = f [v'theta'] / [theta]_p + ...
    #
    ep2_cart = fhat*vertEddy # [1/s*m.hPa/s] = [m.hPa/s2]
    if w is not None:
        w = GetAnomaly(w) # w = w' [hPa/s]
        if wave<0:
            w = np.nanmean(w*u,axis=-1) # w = bar(u'w') [m.hPa/s2]
        else:
            w = GetWaves(u,w,wave=wave) # w = bar(u'w') [m.hPa/s2]
        ep2_cart = ep2_cart - w # [m.hPa/s2]
    #
    #
    # We now have to make sure we get the geometric terms right
    # With our definition,
    #  div1 = 1/(a.cosphi)*d/dphi[a*cosphi*ep1_cart*cosphi],
    #    where a*cosphi comes from using cartesian, and cosphi from the derivative
    # With some algebra, we get
    #  div1 = cosphi d/d phi[ep1_cart] - 2 sinphi*ep1_cart
    div1 = coslat*gradient(ep1_cart,1,1,dphi,edge_order=2)[-1] - 2*sinlat*ep1_cart
    # Now, we want acceleration, which is div(F)/a.cosphi [m/s2]
    div1 = R*div1 # [m/s2]
    #
    # Similarly, we want acceleration = 1/a.coshpi*a.cosphi*d/dp[ep2_cart] [m/s2]
    div2 = gradient(ep2_cart,1,dp,1,edge_order=2)[1] # [m/s2]
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
    from numpy import trapz,cos,prod,reshape,newaxis,pi,where
    #get data into the correct shape
    tmp = AxRoll(data,axis)
    shpe= tmp.shape
    tmp = reshape(tmp,(shpe[0],prod(shpe[1:])))
    #cosine weighting
    J = where((lat>=lim)*(lat<=mx))[0]
    coslat = cos(lat*pi/180.)**cosp
    coswgt = trapz(coslat[J],lat[J])
    tmp = trapz(tmp[J,:]*coslat[J][:,newaxis],lat[J],axis=0)/coswgt
    integ = reshape(tmp,shpe[1:])
    return integ

##############################################################################################
def ComputeN2(pres,Tz,H=7.e3,Rd=287.04,cp=1004):
    ''' Compute the Brunt-Vaisala frequency from zonal mean temperature
         N2 = -Rd*p/(H**2.) * (dTdp - Rd*Tz/(p*cp))
         this is equivalent to
         N2 = g/\theta d\theta/dz, with p = p0 exp(-z/H)

        INPUTS:
            pres  - pressure [hPa]
            Tz    - zonal mean temperature [K], dim pres x lat
            H     - scale height [m]
            Rd    - specific gas constant for dry air
            cp    - specific heat of air at constant pressure
        OUTPUTS:
            N2  - Brunt-Vaisala frequency, [1/s2], dim pres x lat
    '''
    from numpy import newaxis,gradient
    dp   = gradient(pres)[:,newaxis]*100.
    dTdp = gradient(Tz,dp,1,edge_order=2)[0]
    p = pres[:,newaxis]*100. # [Pa]
    N2 = -Rd*p/(H**2.) * (dTdp - Rd*Tz/(p*cp))
    return N2

def ComputeMeridionalPVGrad(lat, pres, uz, Tz, Rd=287.04, cp=1004, a0=6.371e6):
    '''Compute the meridional gradient of potential vorticity.
        This quantity has three terms,
        q_\phi = A - B + C, where
                A = 2*Omega*cos\phi
                B = \partial_\phi[\partial_\phi(ucos\phi)/acos\phi]
                C = af^2/Rd*\partial_p(p\theta\partial_pu/(T\partial_p\theta))

        INPUTS:
            lat  - latitude [degrees]
            pres - pressure [hPa]
            uz   - zonal mean zonal wind [m/s], dim pres x lat OR N x pres x lat
            Tz   - zonal mean temperature [K], dim pres x lat OR N x pres x lat
        OUTPUTS:
            q_phi - meridional gradient of potential vorticity [1/s], dim pres x lat OR N x pres x lat
    '''
    from numpy import pi,cos,sin,newaxis,gradient,deg2rad
    # some constants
    Omega = 2*pi/(86400.) # [1/s]
    p0    = 1e5 #[Pa]

    ## make sure we have the dimesions as expected
    if uz.shape != Tz.shape:
        raise ValueError('UZ AND TZ DO NOT HAVE THE SAME SHAPE')
    elif len(uz.shape) > 3:
        raise ValueError('TOO MANY DIMENSIONS IN UZ AND TZ')
    def FlexiGradPhi(data,dphi):
        if len(data.shape) == 3:
            grad = gradient(data,1,1,dphi,edge_order=2)[2]
        else:
            grad = gradient(data,1,dphi,edge_order=2)[1]
        return grad
    def FlexiGradP(data,dp):
        if len(data.shape) == 3:
            grad = gradient(data,1,dp,1,edge_order=2)[1]
        else:
            grad = gradient(data,dp,1,edge_order=2)[0]
        return grad

    ## convert to Pa
    p = pres[:]*100
    if len(uz.shape) == 3:
        dp = gradient(p)[newaxis,:,newaxis]
        p  = p[newaxis,:,newaxis]
    else:
        dp = gradient(p)[:,newaxis]
        p = p[:,newaxis]
    ## convert to radians
    latpi = deg2rad(lat)
    if len(uz.shape) == 3:
        dphi = gradient(latpi)[newaxis,newaxis,:]
        latpi= latpi[newaxis,newaxis,:]
    else:
        dphi = gradient(latpi)[newaxis,:]
        latpi = latpi[newaxis,:]

    #
    ## first term A
    A = 2*Omega*cos(latpi)

    #
    ## second term B
    dudphi = FlexiGradPhi(uz*cos(latpi),dphi)
    # dudphi = gradient(uz*cos(latpi),1,dphi,edge_order=2)[1]
    B = dudphi/cos(latpi)/a0
    B = FlexiGradPhi(B,dphi)
    # B = gradient(B,1,dphi,edge_order=2)[1]

    #
    ## third term C
    f = 2*Omega*sin(latpi)

    dudp = FlexiGradP(uz,dp)
    # dudp = gradient(uz,dp,1,edge_order=2)[0]

    kappa = Rd/cp
    pp0   = (p0/p)**kappa
    theta = Tz*pp0
    theta_p = FlexiGradP(theta,dp)
    # theta_p = gradient(theta,dp,1,edge_order=2)[0]

    C = p*theta*dudp/(Tz*theta_p)
    C = FlexiGradP(C,dp)
    # C = gradient(C,dp,1,edge_order=2)[0]
    C = a0*f*f*C/Rd

    return A-B+C


def ComputeRefractiveIndex(lat,pres,uz,Tz,k,N2const=None):
    '''
        Refractive index as in Simpson et al (2009) doi 10.1175/2008JAS2758.1 and also Matsuno (1970) doi 10.1175/1520-0469(1970)027<0871:VPOSPW>2.0.CO;2
        Stationary waves are assumed, ie c=0.

        meridonal PV gradient is
        q_\phi = A - B + C, where
                A = 2*Omega*cos\phi
                B = \partial_\phi[\partial_\phi(ucos\phi)/acos\phi]
                C = af^2/Rd*\partial_p(p\theta\partial_pu/(T\partial_p\theta))
        Total refractive index is
        n2 = a^2*[D - E - F], where
                D = q_\phi/(au)
                E = (k/acos\phi)^2
                F = (f/2NH)^2

        Inputs are:
            lat   - latitude [degrees]
            pres  - pressure [hPa]
            uz    - zonal mean zonal wind, dimension pres x lat [m/s]
            Tz    - zonal mean temperature, dimension pres x lat [K]
            k     - zonal wave number [.]
            N2const - if not None, assume N2 = const = N2const [1/s2]
        Outputs are:
            n2  - refractive index, dimension pres x lat [.]
    '''
    from numpy import cos,sin,deg2rad
    # some constants
    Rd    = 287.04 # [J/kg.K = m2/s2.K]
    cp    = 1004 # [J/kg.K = m2/s2.K]
    a0    = 6.371e6 # [m]
    Omega = 2*pi/(24*3600.) # [1/s]
    H     = 7.e3 # [m]

    latpi = deg2rad(lat)

    #
    ## term D
    dqdy = ComputeMeridionalPVGrad(lat,pres,uz,Tz,Rd,cp,a0)
    D = dqdy/(a0*uz)

    #
    ## term E
    latpi = latpi[newaxis,:]
    E = ( k/(a0*cos(latpi)) )**2

    #
    ## term F
    f = 2*Omega*sin(latpi)
    f2 = f*f
    if N2const is None:
        N2 = ComputeN2(pres,Tz,H,Rd,cp)
    else:
        N2 = N2const
    H2 = H*H
    F = f2/(4*N2*H2)

    return a0*a0*(D-E-F)

##############################################################################################
def GetWaves(x,y=[],wave=-1,axis=-1,do_anomaly=False):
	"""Get Fourier mode decomposition of x, or <x*y>, where <.> is zonal mean.

        If y!=[], returns Fourier mode contributions (amplitudes) to co-spectrum zonal mean of x*y. Shape is same as input, except axis which is len(axis)/2+1 due to Fourier symmetry for real signals.

        If y=[] and wave>=0, returns real space contribution of given wave mode. Output has same shape as input.
        If y=[] and wave=-1, returns real space contributions for all waves. Output has additional first dimension corresponding to each wave.

	INPUTS:
		x          - the array to decompose
		y          - second array if wanted
		wave       - which mode to extract. all if <0
		axis       - along which axis of x (and y) to decompose
		do_anomaly - decompose from anomalies or full data
	OUTPUTS:
		xym        - data in Fourier space
	"""
	from numpy import fft,squeeze,real,zeros,zeros_like
	initShape = x.shape
	x = AxRoll(x,axis)
	# compute anomalies
	if do_anomaly:
            x = GetAnomaly(x,0)
        if len(y) > 0:
            y = AxRoll(y,axis)
            if do_anomaly:
                y = GetAnomaly(y,0)
    # Fourier decompose
	x = fft.fft(x,axis=0)
	nmodes = x.shape[0]/2+1
	if wave < 0:
		if len(y) > 0:
			xym = zeros((nmodes,)+x.shape[1:])
		else:
			xym = zeros((nmodes,)+initShape)
	if len(y) > 0:
            y = fft.fft(y,axis=0)
            # Take out the waves
            nl  = x.shape[0]**2
            xyf  = real(x*y.conj())/nl
            # due to symmetric spectrum, there's a factor of 2, but not for wave-0
            mask = zeros_like(xyf)
            if wave < 0:
            	for m in range(xym.shape[0]):
            		mask[m,:] = 1
            		mask[-m,:]= 1
            		xym[m,:] = sum(xyf*mask,axis=0)
            		mask[:] = 0
            	xym = AxRoll(xym,axis,'i')
            else:
            	xym = xyf[wave,:]
            	if wave >= 0:
                	xym = xym + xyf[-wave,:]
	else:
            mask = zeros_like(x)
            if wave >= 0:
                mask[wave,:] = 1
                mask[-wave,:]= 1 # symmetric spectrum for real signals
                xym = real(fft.ifft(x*mask,axis=0))
                xym = AxRoll(xym,axis,'i')
            else:
                for m in range(xym.shape[0]):
                    mask[m,:] = 1
                    mask[-m,:]= 1 # symmetric spectrum for real signals
                    fourTmp = real(fft.ifft(x*mask,axis=0))
                    xym[m,:] = AxRoll(fourTmp,axis,'i')
                    mask[:] = 0
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
    a0    = 6.371e6
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
    N = np.where( (lat>= minLat)*(lat<= maxLat) )[0]
    #  make sure first index is closer to pole
    if lat[N[0]] < lat[N[-1]]:
        latN = [N[0],N[-1]]
    else:
        latN = [N[-1],N[0]]
    S = np.where( (lat<=-minLat)*(lat>=-maxLat) )[0]
    #  make sure first index is closer to pole
    if lat[S[0]] > lat[S[-1]]:
        latS = [S[0],S[-1]]
    else:
        latS = [S[-1],S[0]]
    # pressure levels
    if pres is None:
        temp = tempIn[:,np.newaxis,:]
        K = [0]
    else:
        temp = tempIn
        K = np.where(pres >= minPres)[0]
        if len(K) == 1:
            temp = tempIn[:,K,:]
            K = [0]
            pres = None
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
            tmp = np.squeeze(T[H])
        else:
            tmp = np.trapz(T[H],pres[K],axis=1)
            tmp = tmp/np.trapz(np.ones_like(pres[K]),pres[K])
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
    return cb

#######################################################
def Convert2Days(time,units,calendar):
    """
    """
    import netCDF4 as nc
    import netcdftime as nct
    date = nc.num2date(time,units,calendar)
    unitArray = units.split()
    dayUnits = units.replace(unitArray[0],'days')
    t = nct.utime(dayUnits,calendar=calendar)
    return t.date2num(date)
