#!/usr/bin/python
# Filename: climate.py
#
# Code by Martin Jucker, distributed under an GPLv3 License
#
# This file provides helper functions that can be useful as pre-viz-processing of files and data

############################################################################################
#
from __future__ import print_function
import numpy as np
# from numba import jit

## helper function: Get actual width and height of axes
def GetAxSize(fig,ax,dpi=False):
	"""get width and height of a given axis.
	   output is in inches if dpi=False, in dpi if dpi=True
	"""
	bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
	width, height = bbox.width, bbox.height
	if dpi:
		width *= fig.dpi
		height *= fig.dpi
	return width, height

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
		file	    file name, relative path from wkdir
		climatType  'daily', 'monthly', 'annual', 'DJF', 'JJA', or any
					combination of months according to two-letter code
					Ja Fe Ma Ap My Jn Jl Au Se Oc No De
		wkdir	    working directory, in which 'file' must be, and to which the output
					is written
		timeDim	    name of the time dimension in the netcdf file
		cal	    calendar, if other than within the netcdf file
		Outputs:
		outFile	    name of the output file created
		writes outputfile with name depending on input file name and climatType
	"""

	# need to read netCDF and of course do some math
	import netCDF4 as nc
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

	time	= ncFile.variables[timeDim][:]
	numTimeSteps = len(time)
	timeVar = ncFile.variables[timeDim]
	# check the time units
	timeUnits = timeVar.units
	chck = CheckAny(timeUnits,('seconds','days','months'))
	if not chck:
		print('Cannot understand units of time, which is: '+timeUnits)
		newUnits = raw_input('Please provide units [seconds,days,months] ')
		if newUnits not in ["seconds","days","months"]:
			raise ValueError('units must be seconds, days, or months')
		unitSplit = timeUnits.split()
		unitSplit[0] = newUnits
		timeUnits = ' '.join(unitSplit)
	timeStep = np.diff(timeVar).mean()
	print('The time dimension is in units of',timeUnits,', with a mean time step of',timeStep,'days')
	# check the calendar type
	getCal = False
	if cal:
		timeCal = cal
	else:
		try:
			timeCal = str(timeVar.calendar)
			if not CheckAny(timeCal,calendar_types):
				print('Cannot understand the calendar type, which is: '+timeCal)
				timeCal = raw_input('Please provide a calendar type from the list '+str(calendar_types)+' ')
				timeVar.calendar = timeCal
		except:
			timeCal = raw_input('Please provide a calendar type from the list '+str(calendar_types)+' ')
	if timeCal not in calendar_types:
		raise ValueError('calender must be in '+str(calendar_types))
	else:
		print('Calendar type '+timeCal)
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

	print('Averaging variables:')
	for var in ncFile.variables:
		varShape = np.shape(ncFile.variables[var])
		if len(varShape) == 0: continue
		if varShape[0] == numTimeSteps and len(varShape) >= 2:
			print('			    ',var)
			tmpVar = ncFile.variables[var][:]
			if climType != 'daily' and climType != 'monthly':
				outVar = outFile.createVariable(var,str(ncFile.variables[var].dtype),ncFile.variables[var].dimensions[1:])
				tmpAvg = tmpVar[climTimeVar>0,:].mean(axis=0)
			else:
				outVar = outFile.createVariable(var,str(ncFile.variables[var].dtype),ncFile.variables[var].dimensions	 )
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
	print('DONE, wrote file',outFileName)
	return outFileName

##############################################################################################
# get the saturation mixing ration according to Clausius-Clapeyron

# helper function: re-arrange array dimensions
def AxRoll(x,ax,invert=False):
	"""Re-arrange array x so that axis 'ax' is first dimension.
		Undo this if invert=True
	"""
	if ax < 0:
		n = len(x.shape) + ax
	else:
		n = ax
	#
	if invert is False:
		y = np.rollaxis(x,n,0)
	else:
		y = np.rollaxis(x,0,n+1)
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
	#some constants we need
	from .constants import Rd,Rv,ESO,HLV,Tfreeze

	# make sure we are operating along the pressure axis
	T = AxRoll(T,pDim)

	# pressure is assumed in hPa: convert to Pa
	p = p*100
	# compute saturation pressure
	esat = ES0*np.exp(HLV*(1./Tfreeze - 1./T)/Rv)
	qsat = np.zeros_like(esat)
	# finally, compute saturation mixing ratio from pressure
	for k in range(len(p)):
		qsat[k,:] = Rd/Rv*esat[k,:]/(p[k]-esat[k,:])
	return AxRoll(qsat,pDim,invert=True)


##############################################################################################
def ComputeRelativeHumidity(inFile, pDim, outFile='none', temp='temp', sphum='sphum', pfull='pfull'):
	"""Computes relative humidity from temperature and specific humidity.

		File inFile is assumed to contain both temperature and specific humidity.
		Relative humidity is either output of the function, or written to the file outFile.

		Inputs:
			inFile	  Name of the file (full path)
						containing temperature and moisture
			pDim	  Index of pressure dimension within temperature array
			outFile	  Name of the output file containing specific humidity.
						No output file is created if outFile='none'
			temp	  Name of the temperature variable inside inFile
			sphum	  Name of specific humidity variable inside inFile
			pfull	  Name of full level pressure [hPa] inside inFile
	"""

	import netCDF4 as nc
	# relative humidity is then q/qsat*100[->%]

	# read input file
	inFile = nc.Dataset(inFile, 'r')
	t = inFile.variables[temp][:]
	q = inFile.variables[sphum][:]
	p = inFile.variables[pfull][:]

	# compute saturation mixing ratio
	qsat = ComputeSaturationMixingRatio(t, p, pDim)

	#write output file
	if outFile != 'none':
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
			data	    - filename of input file or dictionary with temp,vcomp,lat,pfull
			outFileName - filename of output file, 'none', or 'same'
			temp	    - name of temperature field in inFile
			vcomp	    - name of meridional velocity field in inFile
			lat	    - name of latitude in inFile
			pfull	    - name of pressure in inFile [hPa]
			time	    - name of time field in inFile. Only needed if outFile used
			p0	    - pressure basis to compute potential temperature [hPa]
		OUTPUTS:
			psi	    - stream function, as a function of time
			psis	    - residual stream function, as a function of time
	"""
	import netCDF4 as nc
	from scipy.integrate import cumtrapz
	import os

	# some constants
	from .constants import kappa,a0,g

	if isinstance(data,str):
		# check if file exists
		if not os.path.isfile(data):
			raise IOError('File '+data+' does not exist')
		# read input file
		print('Reading data')
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
	p  = p *100 # [Pa]
	p0 = p0*100 # [Pa]
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
	if outFileName != 'none':
		print('Writing file  '+outFileName)
		if outFileName != 'same':
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
		print('Done writing file '+outFileName)
		if outFileName != 'same':
			inFile.close()
	return psi,psis

##############################################################################################
def ComputePsiXr(v, t, lon='lon', lat='lat', pres='level', time='time', ref='mean', p0=1e3):
	"""Computes the residual stream function \Psi* (as a function of time).

		INPUTS:
			v	    - meridional wind, xr.DataArray
			t	    - temperature, xr.DataArray
			lon	    - name of longitude in t
			lat	    - name of latitude in t
			pres	    - name of pressure in t [hPa]
			time	    - name of time field in t
			ref	    - how to treat dTheta/dp:
				       - 'rolling-X' : centered rolling mean over X days
				       - 'mean'	     : full time mean
			p0	    - pressure basis to compute potential temperature [hPa]
		OUTPUTS:
			psi	    - stream function, as a function of time
			psis	    - residual stream function, as a function of time
	"""
	from scipy.integrate import cumtrapz
	from numpy import cos,deg2rad
	# some constants
	from .constants import kappa,a0,g
	#
	## compute psi
	v_bar,t_bar = ComputeVertEddyXr(v,t,pres,p0,lon,time,ref) # t_bar = bar(v'Th'/(dTh_bar/dp))

	# Eulerian streamfunction
	pdim = v_bar.get_axis_num(pres)
	psi = v_bar.reduce(cumtrapz,x=v_bar[pres],axis=pdim,initial=0) # [m.hPa/s]


	## compute psi* = psi - bar(v'Th'/(dTh_bar/dp))
	psis = psi - t_bar
	coslat = np.cos(np.deg2rad(t[lat]))
	psi = 2*np.pi*a0/g*psi *coslat*100 #[kg/s]
	psis= 2*np.pi*a0/g*psis*coslat*100 #[kg/s]

	## give the DataArrays their names
	psi.name = 'psi'
	psis.name= 'psis'

	return psi,psis


##############################################################################################
## helper functions
def update_progress(progress,barLength=10,info=None):
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
		status = "\r" #"\r\n"
	#status = "Done...\r\n"
	block = int(round(barLength*progress))
	if info is not None:
		text = '\r'+info+': '
	else:
		text = '\r'
	if progress == 1:
		if info is not None:
			text = "\r{0}	{1}	{2}".format(" "*(len(info)+1)," "*barLength,status)
		else:
			text = "\r   {0}     {1}".format(" "*barLength,status)
	else:
		text += "[{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), int(progress*100), status)
	sys.stdout.write(text)
	sys.stdout.flush()
#
def ComputeVertEddy(v,t,p,p0=1e3,wave=0):
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
	#
	# some constants
	from .constants import kappa
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
	dthdp = np.gradient(t_bar,edge_order=2)[1]/dp # dthdp = d(theta_bar)/dp
	dthdp[dthdp==0] = np.NaN
	# time mean of d(theta_bar)/dp
	dthdp = np.nanmean(dthdp,axis=0)[np.newaxis,:]
	# now get wave component
	#if isinstance(wave,list):
	#	t = np.sum(GetWaves(v,t,wave=-1,do_anomaly=True)[:,:,:,wave],axis=-1)
	if wave == 0:
		v = GetAnomaly(v) # v = v'
		t = GetAnomaly(t) # t = t'
		t = np.nanmean(v*t,axis=-1) # t = bar(v'Th')
	else:
		t = GetWaves(v,t,wave=wave,do_anomaly=True) # t = bar(v'Th'_{k=wave})
		if wave < 0:
			dthdp = np.expand_dims(dthdp,-1)
	t_bar = t/dthdp # t_bar = bar(v'Th')/(dTh_bar/dp)
	#
	return v_bar,t_bar

def ComputeVertEddyXr(v,t,p='level',p0=1e3,lon='lon',time='time',ref='mean',wave=0):
	""" Computes the vertical eddy components of the residual circulation,
		bar(v'Theta'/Theta_p).
		Output units are [v_bar] = [v], [t_bar] = [v*p]

		INPUTS:
			v    - meridional wind, xr.DataArray
			t    - temperature, xr.DataArray
			p    - name of pressure
			p0   - reference pressure for potential temperature
			lon  - name of longitude
			time - name of time field in t
			ref  - how to treat dTheta/dp:
			       - 'rolling-X' : centered rolling mean over X days
			       - 'mean'	     : full time mean
                               - 'instant'   : no time operation
			wave - wave number: if == 0, return total. else passed to GetWavesXr()
		OUPUTS:
			v_bar - zonal mean meridional wind [v]
			t_bar - zonal mean vertical eddy component <v'Theta'/Theta_p> [v*p]
	"""
	#
	# some constants
	from .constants import kappa
	#
	# pressure quantitites
	pp0 = (p0/t[p])**kappa
	# convert to potential temperature
	t = t*pp0 # t = theta
	# zonal means
	v_bar = v.mean(lon)
	t_bar = t.mean(lon) # t_bar = theta_bar
	# prepare pressure derivative
	dthdp = t_bar.differentiate(p,edge_order=2) # dthdp = d(theta_bar)/dp
	dthdp = dthdp.where(dthdp != 0)
	# time mean of d(theta_bar)/dp
	if time in dthdp.dims:
		if 'rolling' in ref:
			r = int(ref.split('-')[-1])
			dthdp = dthdp.rolling(dim={time:r},min_periods=1,center=True).mean()
		elif ref == 'mean':
			dthdp = dthdp.mean(time)
		elif ref == 'instant':
			dthdp = dthdp
	# now get wave component
	if isinstance(wave,list):
		vpTp = GetWavesXr(v,t,dim=lon,wave=-1).sel(k=wave).sum('k')
	elif wave == 0:
		vpTp = (v - v_bar)*(t - t_bar)
		vpTp = vpTp.mean(lon)  # vpTp = bar(v'Th')
	else:
		vpTp = GetWavesXr(v,t,dim=lon,wave=wave) # vpTp = bar(v'Th'_{k=wave})
	t_bar = vpTp/dthdp # t_bar = bar(v'Th')/(dTh_bar/dp)
	#
	return v_bar,t_bar

##############################################################################################
def eof(X,n=-1,detrend='constant',eof_in=None):
	"""Principal Component Analysis / Empirical Orthogonal Functions / SVD

		Uses Singular Value Decomposition to find the dominant modes of variability.
		The field X can be reconstructed with Y = dot(EOF,PC) + X.mean(axis=time)

		INPUTS:
			X	-- Field, shape (time x space).
			n	-- Number of modes to extract. All modes if n < 0
			detrend -- detrend with global mean ('constant')
						  or linear trend ('linear')
		    eof_in  -- If not None, compute PC by projecting eof onto X.
		OUTPUTS:
			EOF - Spatial modes of variability
			PC  - Temporal evolution of EOFs - only output if eof_in is not None
			E   - Explained value of variability
			u   - spatial modes
			s   - variances
			v   - temporal modes
	"""
	is_xr = False
	try:
		import xarray as xr
		if isinstance(X,xr.DataArray):
			is_xr = True
			dims = []
			for dim in X.dims[1:]:
				dims.append(X[dim])
			tdim = [X[X.dims[0]]]
			X = X.values
	except:
		pass
	import scipy.signal as sg
	# make sure we have a matrix time x space
	shpe = X.shape
	if len(shpe) > 2:
		X = X.reshape([shpe[0],np.prod(shpe[1:])])
		if eof_in is not None:
			if len(eof_in.shape) > 2:
				eof_in = eof_in.reshape([np.prod(eof_in.shape[:-1]),eof_in.shape[-1]])
			else:
				eof_in = eof_in.reshape([np.prod(eof_in.shape),1])
	# take out the time mean or trend
	X = sg.detrend(X.transpose(),type=detrend)
	if eof_in is not None:
		if eof_in.shape[-1] == X.shape[0]:
			PC =  np.matmul(eof_in, X)
			eof_norm = np.dot(eof_in.transpose(),eof_in)
			return np.dot(PC,np.linalg.inv(eof_norm))
		else:
			PC = np.matmul(eof_in.transpose(), X)
			eof_norm = np.dot(eof_in.transpose(),eof_in)
			return np.dot(PC.transpose(),np.linalg.inv(eof_norm)).transpose()
		# return sg.detrend(PC,type='constant')
	# perform SVD - v is actually V.H in X = U*S*V.H
	u,s,v = np.linalg.svd(X, full_matrices=False)
	# now, u contains the spatial, and v the temporal structures
	# s contains the variances, with the same units as the input X
	# u.shape = (space, modes(space)), v.shape = (modes(space), time)

	# get the first n modes, in physical units
	#  we can either project the data onto the principal component, X*V
	#  or multiply u*s. This is the same, as U*S*V.H*V = U*S
	if n < 0:
		n = s.shape[0]
	EOF = np.dot(u[:,:n],np.diag(s)[:n,:n])
	# time evolution is in v
	PC  = v[:n,:]
	# EOF wants \lambda = the squares of the eigenvalues,
	#  but SVD yields \gamma = \sqrt{\lambda}
	s2 = s*s
	E   = s2[:n]/sum(s2)
	# now we need to make sure we get everything into the correct shape again
	u = u[:,:n]
	s = s[:n]
	v = v.transpose()[:,:n]
	if len(shpe) > 2:
		# replace time dimension with modes at the end of the array
		newshape = list(shpe[1:])+[n]
		EOF = EOF.reshape(newshape)
		u   = u	 .reshape(newshape)
	if is_xr: # return xarray dataarrays
		mode = [('n',np.arange(1,n+1))]
		EOF = xr.DataArray(EOF,coords=dims+mode,name='EOF')
		PC  = xr.DataArray(PC,coords=tdim+mode,name='PC')
		E   = xr.DataArray(E,coords=mode,name='E')
	return EOF,PC,E,u,s,v


##############################################################################################
def ComputeAnnularMode(lat, pres, data, choice='z', hemi='infer', detrend='constant', eof_in=None, pc_in=None, eof_out=False, pc_out=False):
	"""Compute annular mode as in Gerber et al, GRL 2008.
		This is basically the first PC, but normalized to unit variance and zero mean.
		To conform to Gerber et al (2008), `data` should be anomalous height or zonal wind
		 with respect to 30-day smoothed day of year climatology.

		INPUTS:
			lat    - latitude
			pres   - pressure
			data   - variable to compute EOF from. This is typically
						geopotential or zonal wind.
						Size time x pres x lat (ie zonal mean)
			choice - not essential, but used for sign convention.
						If 'z', the sign is determined based on 70-80N/S.
						Otherwise, 50-60N/S is used.
			hemi   - hemisphere to consider
						'infer' - if mean(lat)>=0 -> NH, else SH
						'SH' or 'NH'
			detrend- detrend method for computing EOFs:
						'linear' -> remove linear trend
						'constant' -> remove total time mean
			eof_in - if None, compute EOF1 as usual.
					 if the EOF1 is already known, use this instead of
				    computing it again.
			pc_in  - if None, standardize PC1 to its own mean and std deviation
				 else, use pc_in mean and std deviation to standardize.
			eof_out- whether or not to pass the first EOF as output [False].
			pc_out - whether or not to pass the first PC as output [False].
		OUTPUT:
			AM     - The annular mode, size time x pres
			EOF    - The first EOF (if eof_out is True), size pres x lat
			PC     - The first PC (if pc_out is True). size time x pres
	"""
	#
	AM = np.full((data.shape[0],data.shape[1]),np.nan)
	if pc_out:
		pco = np.full(AM.shape,np.nan)
	# guess the hemisphere
	if hemi == 'infer':
		if np.mean(lat) >= 0:
			sgn = 1.
		else:
			sgn = -1.
	elif hemi == 'SH':
		sgn = -1.
	elif hemi == 'NH':
		sgn = 1.
	j_tmp = np.where(sgn*lat > 20)[0]
	if eof_out:
		eofo = np.full((data.shape[1],len(j_tmp)),np.nan)
	coslat = np.cos(np.deg2rad(lat))
	negCos = (coslat < 0.)
	coslat[negCos] = 0.
	# weighting as in Gerber et al GRL 2008
	sqrtcoslat = np.sqrt(coslat[j_tmp])
	# try to get the sign right
	# first possibility
	if choice == 'z':
		minj = min(sgn*70,sgn*80)
		maxj = max(sgn*80,sgn*70)
		sig = -1
	else:
		minj = min(sgn*50,sgn*60)
		maxj = max(sgn*60,sgn*50)
		sig = 1
	jj = (lat[j_tmp] > minj)*(lat[j_tmp] < maxj)
	# second possibility
	#jj = abs(lat[j_tmp]-80).argmin()
	#sig = -1
	if isinstance(pres,(int,float)):
		data = np.reshape(data,(data.shape[0],1,data.shape[1]))
		pres = [pres]
	for k in range(len(pres)):
		# remove global mean
		globZ = GlobalAvg(lat,data[:,k,:],axis=-1,lim=lat[j_tmp[0]],mx=lat[j_tmp[-1]])
		var = data[:,k,:] - globZ[:,np.newaxis]
		# area weighting: EOFs are ~variance, thus take sqrt(cos)
		var = var[:,j_tmp]*sqrtcoslat[np.newaxis,:]
		varNan = np.isnan(var)
		if np.sum(np.reshape(varNan,(np.size(varNan),)))==0:
			if eof_in is None:
				eof1,pc1,E,u,s,v = eof(var,n=1,detrend=detrend)
			else:
				pc1 = eof(var,n=1,detrend=detrend,eof_in=np.expand_dims(eof_in[k,:],-1))
				eof1 = eof_in[k,:]
			# force the sign of PC
			pc1  = pc1*sig*np.sign(eof1[jj].mean())
			if eof_out:
				eofo[k,:] = np.squeeze(eof1)
			if pc_out:
				pco[:,k] = pc1
			# force unit variance and zero mean
			if pc_in is None:
				AM[:,k] = (pc1-pc1.mean())/np.std(pc1)
			else:
				AM[:,k] = (pc1-pc_in.mean())/np.std(pc_in)
	if eof_out and pc_out:
		return AM,eofo,pco
	elif eof_out:
		return AM,eofo
	elif pc_out:
		return AM,pco
	else:
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
			vstar	    - residual meridional wind, as a function of time
	"""
	import netCDF4 as nc

	a0    = 6371000
	g     = 9.81

	# read input file
	if isinstance(data,str):
		print('Reading data')
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
	vstar = v_bar - np.gradient(t_bar,edge_order=2)[1]/dp

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
		vt_bar[w,:] = np.gradient(vt_bar[w,:],edge_order=2)[-1]/dphi
	# compute zonal mean upwelling
	w_bar = np.nanmean(data[omega],axis=-1)
	# put it all together
	if len(wave)==1:
		return w_bar + np.squeeze(R*vt_bar)
	else:
		return w_bar + R*vt_bar

##############################################################################################
def ComputeWstarXr(omega, temp, vcomp, pres='level', lon='lon', lat='lat', time='time', ref='mean', p0=1e3, is_Pa='omega'):
	"""Computes the residual upwelling w*. omega, temp, vcomp are xarray.DataArrays.

		Output units are the same as the units of omega, and the pressure coordinate is expected in hPa, latitude in degrees.

		INPUTS:
			omega - pressure velocity. xarray.DataArray
			temp  - temperature. xarray.DataArray
			vcomp - meridional velocity. xarray.DataArray
			pfull - name of pressure coordinate.
			lon   - name of longitude coordinate
			lat   - name of latitude coordinate
			time  - name of time coordinate
			ref   - how to treat dTheta/dp:
				- 'rolling-X' : centered rolling mean over X days
				- 'mean'      : full time mean
			p0    - pressure basis to compute potential temperature [hPa]
			is_Pa - correct for pressure units in variables:
				- None: omega, p0 and pres are all in hPa or all in Pa
				- 'omega': omega is in Pa/s but pres and p0 in hPa
				- 'pres' : omega is in hPa/s, p0 in hPa, but pres in Pa
		OUTPUTS:
			residual pressure velocity, same units as omega
	"""
	import numpy as np

	a0    = 6371000.

	# spherical geometry
	coslat = np.cos(np.deg2rad(omega[lat]))
	R = a0*coslat
	R = 1./R
	# correct for units: hPa<->Pa
	if is_Pa is not None:
		if is_Pa.lower() == 'omega':
			R = R*100
		elif is_Pa.lower() == 'pres':
			R  = R*0.01
			p0 = p0*100
        # correct for units: degrees<->radians
	R = R*180/np.pi
	# compute thickness weighted meridional heat flux
	_,vt_bar = ComputeVertEddyXr(vcomp, temp, pres, p0, lon, time, ref)
	# get the meridional derivative
	vt_bar = (coslat*vt_bar).differentiate(lat)
	# compute zonal mean upwelling
	w_bar = omega.mean(lon)
	# put it all together
	return w_bar + R*vt_bar

##############################################################################################
def ComputeEPfluxDiv(lat,pres,u,v,t,w=None,do_ubar=False,wave=0):
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
	  wave - only include this wave number. total if == 0, all waves if <0, single wave if >0, sum over waves if a list. optional
	OUTPUTS:
	  ep1  - meridional EP-flux component, scaled to plot in cartesian [m2/s2]
	  ep2  - vertical   EP-flux component, scaled to plot in cartesian [hPa*m/s2]
	  div1 - horizontal EP-flux divergence, divided by acos\phi [m/s/d]
	  div2 - horizontal EP-flux divergence , divided by acos\phi [m/s/d]
	"""
	# some constants
	from .constants import Rd,cp,kappa,p0,Omega,a0
	# geometry
	pilat = lat*np.pi/180
	dphi  = np.gradient(pilat)[np.newaxis,np.newaxis,:]
	coslat= np.cos(pilat)[np.newaxis,np.newaxis,:]
	sinlat= np.sin(pilat)[np.newaxis,np.newaxis,:]
	R     = 1./(a0*coslat)
	f     = 2*Omega*sinlat
	pp0  = (p0/pres[np.newaxis,:,np.newaxis])**kappa
	dp    = np.gradient(pres)[np.newaxis,:,np.newaxis]
	if wave < 0:
		dp   = np.expand_dims(dp  ,-1)
	#
	# absolute vorticity
	if do_ubar:
		ubar = np.nanmean(u,axis=-1)
		fhat = R*np.gradient(ubar*coslat,edge_order=2)[-1]/dphi
	else:
		fhat = 0.
	fhat = f - fhat # [1/s]
	#
	# add wavenumber dimension if needed
	if wave < 0:
		fhat = np.expand_dims(fhat,-1)
		coslat = np.expand_dims(coslat,-1)
		sinlat = np.expand_dims(sinlat,-1)
		dphi = np.expand_dims(dphi,-1)
		R = np.expand_dims(R,-1)
	#
	## compute thickness weighted heat flux [m.hPa/s]
	vbar,vertEddy = ComputeVertEddy(v,t,pres,p0,wave) # vertEddy = bar(v'Th'/(dTh_bar/dp))
	#
	## get zonal anomalies
	u = GetAnomaly(u)
	v = GetAnomaly(v)
	#if isinstance(wave,list):
	#	upvp = np.sum(GetWaves(u,v,wave=-1)[:,:,:,wave],-1)
	if wave == 0:
		upvp = np.nanmean(u*v,axis=-1)
	else:
		upvp = GetWaves(u,v,wave=wave)
	#
	## compute the horizontal component
	if do_ubar:
		grad_u= np.gradient(ubar,edge_order=2)[1]
		if wave < 0:
			grad_u = np.expand_dims(grad_u,-1)
		shear = grad_u/dp # [m/s.hPa]
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
		#if isinstance(wave,list):
		#	w = sum(GetWaves(u,w,wave=wave)[:,:,:,wave],-1)
		if wave == 0:
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
	if wave < 0:
		div1 = coslat*np.gradient(ep1_cart,edge_order=2)[-2]/dphi - 2*sinlat*ep1_cart
	else:
		div1 = coslat*np.gradient(ep1_cart,edge_order=2)[-1]/dphi - 2*sinlat*ep1_cart
	# Now, we want acceleration, which is div(F)/a.cosphi [m/s2]
	div1 = R*div1 # [m/s2]
	#
	# Similarly, we want acceleration = 1/a.coshpi*a.cosphi*d/dp[ep2_cart] [m/s2]
	div2 = np.gradient(ep2_cart,edge_order=2)[1]/dp # [m/s2]
	#
	# convert to m/s/day
	div1 = div1*86400
	div2 = div2*86400
	#
	return ep1_cart,ep2_cart,div1,div2

##############################################################################################
def ComputeEPfluxDivXr(u,v,t,lon='infer',lat='infer',pres='infer',time='time',ref='mean',w=None,do_ubar=False,wave=0):
	""" Compute the EP-flux vectors and divergence terms.

		The vectors are normalized to be plotted in cartesian (linear)
		coordinates, i.e. do not incluxde the geometric factor a*cos\phi.
		Thus, ep1 is in [m2/s2], and ep2 in [hPa*m/s2].
		The divergence is in units of m/s/day, and therefore represents
		the deceleration of the zonal wind. This is actually the quantity
		1/(acos\phi)*div(F).

	INPUTS:
	  u    - zonal wind, xarray.DataArray [m/s]
	  v    - meridional wind, xarray.DataArray [m/s]
	  t    - temperature, xarray.DataArray [K]
	  lon  - name of longitude
	  lat  - name of latitude
	  pres - name of pressure coordinate [hPa/s]
	  time - name of time coordinate (used for ComputeVertEddy)
	  ref  - method to use for dTheta/dp, used for ComputeVertEddy
	  w    - pressure velocity, if not None, xarray.DataArray [hPa/s]
	  do_ubar - compute shear and vorticity correction?
	  wave - only include this wave number. total if == 0, all waves if <0, sum over waves if a list. optional
	OUTPUTS (all xarray.DataArray):
	  ep1  - meridional EP-flux component, scaled to plot in cartesian [m2/s2]
	  ep2  - vertical   EP-flux component, scaled to plot in cartesian [hPa*m/s2]
	  div1 - horizontal EP-flux divergence, divided by acos\phi [m/s/d]
	  div2 - horizontal EP-flux divergence , divided by acos\phi [m/s/d]
	"""
	# some constants
	from .constants import Rd,cp,kappa,p0,Omega,a0
	# shape
	dim_names = FindCoordNames(u)
	if lon == 'infer':
		lon = dim_names['lon']
	if lat == 'infer':
		lat = dim_names['lat']
	if pres == 'infer':
		pres = dim_names['pres']
	initial_order = u.dims
	# geometry
	coslat = np.cos(np.deg2rad(u[lat]))
	sinlat = np.sin(np.deg2rad(u[lat]))
	R      = 1./(a0*coslat)
	f      = 2*Omega*sinlat
	pp0    = (p0/u[pres])**kappa
	#
	# absolute vorticity
	if do_ubar:
		ubar = u.mean(lon)
		fhat = R*np.rad2deg((ubar*coslat)).differentiate(lat,edge_order=2)
	else:
		fhat = 0.
	fhat = f - fhat # [1/s]
	#
	## compute thickness weighted heat flux [m.hPa/s]
	vbar,vertEddy = ComputeVertEddyXr(v,t,pres,p0,lon,time,ref,wave) # vertEddy = bar(v'Th'/(dTh_bar/dp))
	#
	## get zonal anomalies
	if isinstance(wave,list):
		upvp = GetWavesXr(u,v,dim=lon,wave=-1).sel(k=wave).sum('k')
	elif wave == 0:
		u = u - u.mean(lon)
		v = v - v.mean(lon)
		upvp = (u*v).mean(lon)
	else:
		upvp = GetWavesXr(u,v,dim=lon,wave=wave)
	#
	## compute the horizontal component
	if do_ubar:
		shear = ubar.differentiate(pres,edge_order=2) # [m/s.hPa]
	else:
		shear = 0.
	ep1_cart = -upvp + shear*vertEddy # [m2/s2 + m/s.hPa*m.hPa/s] = [m2/s2]
	#
	## compute vertical component of EP flux.
	## at first, keep it in Cartesian coordinates, ie ep2_cart = f [v'theta'] / [theta]_p + ...
	#
	ep2_cart = fhat*vertEddy # [1/s*m.hPa/s] = [m.hPa/s2]
	if w is not None:
		if isinstance(wave,list):
			w = GetWavesXr(u,w,dim=lon,wave=-1).sel(k=wave).sum('k')
		elif wave == 0:
			w = w - w.mean(lon) # w = w' [hPa/s]
			w = (w*u).mean(lon) # w = bar(u'w') [m.hPa/s2]
		else:
			w = GetWavesXr(u,w,dim=lon,wave=wave)
		ep2_cart = ep2_cart - w # [m.hPa/s2]
	#
	#
	# We now have to make sure we get the geometric terms right
	# With our definition,
	#  div1 = 1/(a.cosphi)*d/dphi[a*cosphi*ep1_cart*cosphi],
	#    where a*cosphi comes from using cartesian, and cosphi from the derivative
	# With some algebra, we get
	#  div1 = cosphi d/d phi[ep1_cart] - 2 sinphi*ep1_cart
	div1 = coslat*(np.rad2deg(ep1_cart).differentiate(lat,edge_order=2)) \
			- 2*sinlat*ep1_cart
	# Now, we want acceleration, which is div(F)/a.cosphi [m/s2]
	div1 = R*div1 # [m/s2]
	#
	# Similarly, we want acceleration = 1/a.coshpi*a.cosphi*d/dp[ep2_cart] [m/s2]
	div2 = ep2_cart.differentiate(pres,edge_order=2) # [m/s2]
	#
	# convert to m/s/day
	div1 = div1*86400
	div2 = div2*86400
	#
	# make sure order is the same as input
	new_order = [d for d in initial_order if d != lon]
	if not isinstance(wave,list) and wave < 0:
		new_order = ['k'] + new_order
	# give the DataArrays their names
	ep1_cart.name = 'ep1'
	ep2_cart.name = 'ep2'
	div1.name = 'div1'
	div2.name = 'div2'
	return ep1_cart.transpose(*new_order),ep2_cart.transpose(*new_order),div1.transpose(*new_order),div2.transpose(*new_order)

##############################################################################################
def ComputeStreamfunction(u,v,lat='lat',lon='lon',use_windspharm=False,lat0=0,lon0=0,method='uv',smooth=None,vw=None):
	'''
		Compute the horizontal streamfunction from u and v. Assumes horizontal non-divergence.
		If windspharm is to be used (and installed), this will be the most accurate. However, due to
		some serious restrictions in terms of domain and grid spacing with spherical harmonics,
		a direct integral method is also implemented. Note that this will give a much noisier streamfunction
		than using windspharm.

		INPUTS:
			u	       : zonal wind. xarray.DataArray.
			v	       : meridional wind. xarray.DataArray.
			lat	       : name of latitude in u,v.
			lon	       : name of longitude on u,v.
			use_windspharm : whether or not to use windspharm. If False, use direct integration method.
			lat0	       : reference latitude for direct integration method.
			lon0	       : reference longitude for direct integration method.
			method	       : integrate over u ('u') or v ('v') or both ('uv'). Only for direct integration method.
			smooth	       : smooth streamfunction after integration via rolling mean. If not None, should be a dictionary
					  containing the roll window for each dimension, e.g. {'lon':5,'lat':3}.
			vw	       : if use_windspharm = True, save time by also sending VectorWind object.
		OUTPUTS:
			psi : streamfunction
			wv  : windspharm.VectorWind object. Only applies if use_windspharm = True
	'''
	import numpy as np
	from xarray import DataArray,concat
	a0 = 6376000
	if use_windspharm:
		if vw is None:
			from windspharm.xarray import VectorWind
			vw = VectorWind(u,v)
		psi = vw.streamfunction()
		# windspharm might invert latitude
		return psi.reindex_like(u)
	else:
		u = u.sortby(lat)
		v = v.sortby(lat)
		cosphi = np.cos(np.deg2rad(u[lat]))
		lons = u[lon].values
		lats = u[lat].values
		nlon = len(lons)
		nlat = len(lats)
		# nlats = len(lats)
		# nlons = len(lons)
		j0 = np.argmin(np.abs(lats-lat0))
		i0 = np.argmin(np.abs(lons-lon0))
		from scipy.integrate import cumtrapz
		#
		## first, fix phi0
		#
		if 'u' in method:
			latind = u.get_axis_num(lat)
			# integrate from phi0 to phi_max
			integrand = a0*u
			if j0 < nlat-1:
				integrand_a = integrand.isel({lat:slice(j0,None)})
				psi_1a = -cumtrapz(integrand_a,x=lats[j0:],axis=latind,initial=0)
				psi1ax = DataArray(psi_1a,coords=integrand_a.coords,name='psi')
			else:
				psi1ax = 0*integrand
			# integrate from phi0 to phi_min
			if j0 > 0:
				j1 = j0+1
				integrand_b = integrand.isel({lat:slice(None,j1)}).isel({lat:slice(None,None,-1)})
				psi_1b = -cumtrapz(integrand_b,x=lats[:j1][::-1],axis=latind,initial=0)
				psi1bx = DataArray(psi_1b,coords=integrand_b.coords,name='psi').isel({lat:slice(1,None)})
			else:
				psi1bx = 0*integrand
			# case 1: either psi1ax = 0 or psi1bx = 0
			if psi1ax.shape == integrand.shape and psi1bx.shape == integrand.shape:
				psi1x = psi1ax + psi1bx
			# case 2: integrated on either side of j0
			else:
				psi1x = concat([psi1ax,psi1bx],dim=lat).sortby(lat)

			# integrate at phi0
			integrand = (a0*cosphi*v).isel({lat:j0})
			lonind = integrand.get_axis_num(lon)
			if i0 < nlon-1:
				integrand_a = integrand.isel({lon:slice(i0,None)})
				psi_2a = cumtrapz(integrand_a,x=lons[i0:],axis=lonind,initial=0)
				psi2ax = DataArray(psi_2a,coords=integrand_a.coords,name='psi')
			else:
				psi2ax = 0*integrand
			if i0 > 0:
				i1 = i0+1
				integrand_b = integrand.isel({lon:slice(None,i1)}).isel({lon:slice(None,None,-1)})
				psi_2b = cumtrapz(integrand_b,x=lons[:i1][::-1],axis=lonind,initial=0)
				psi2bx = DataArray(psi_2b,coords=integrand_b.coords,name='psi').isel({lon:slice(1,None)})
			else:
				psi2bx = 0*integrand
			# case 1: either psi2ax = 0 or psi2bx = 0
			if psi2ax.shape == integrand.shape or psi2bx.shape == integrand.shape:
				psi2x = psi2ax + psi2bx
			# case 2: integrated on either side of i0
			else:
				psi2x = concat([psi2ax,psi2bx],dim=lon).sortby(lon)
			# put everything together
			psix1 = psi2x + psi1x
			if method == 'u':
				psi_out = psix1*np.pi/180
		#
		## then, do the same but fix lambda0
		#
		# integrate from lambda0 to lambda
		if 'v' in method:
			integrand = a0*cosphi*v
			lonind = integrand.get_axis_num(lon)
			if i0 < nlon-1:
				integrand_a = integrand.isel({lon:slice(i0,None)})
				psi_1a = cumtrapz(integrand_a,x=lons[i0:],axis=lonind,initial=0)
				psi1ax	= DataArray(psi_1a,integrand_a.coords,name='psi')
			else:
				psi1ax = 0*integrand
			if i0 > 0:
				i1 = i0+1
				integrand_b = integrand.isel({lon:slice(None,i1)}).isel({lon:slice(None,None,-1)})
				psi_1b = cumtrapz(integrand_b,x=lons[:i1][::-1],axis=lonind,initial=0)
				psi1bx = DataArray(psi_1b,integrand_b.coords,name='psi').isel({lon:slice(1,None)})
			else:
				psi1bx = 0*integrand
			if psi1ax.shape == integrand.shape or psi1bx.shape == integrand.shape:
				psi1x = psi1ax + psi1bx
			else:
				psi1x = concat([psi1ax,psi1bx],dim=lon).sortby(lon)
			# them integrate at lambda0
			integrand = a0*u
			latind = integrand.get_axis_num(lat)
			if j0 < nlat-1:
				integrand_a = integrand.isel({lat:slice(j0,None),lon:i0})
				psi_2a = cumtrapz(integrand_a,x=lats[j0:],axis=latind,initial=0)
				psi2ax = DataArray(psi_2a,coords=integrand_a.coords,name='psi')
			else:
				psi2ax = 0*integrand
			if j0 > 0:
				j1 = j0+1
				integrand_b = integrand.isel({lat:slice(None,j1),lon:i0}).isel({lat:slice(None,None,-1)})
				psi_2b = cumtrapz(integrand_b,x=lats[:j1][::-1],axis=latind,initial=0)
				psi2bx = DataArray(psi_2b,coords=integrand_b.coords,name='psi').isel({lat:slice(1,None)})
			else:
				psi2bx = 0*integrand
			if psi2ax.shape == integrand.shape or psi2bx.shape == integrand.shape:
				psi2x = psi2ax + psi2bx
			else:
				psi2x = concat([psi2ax,psi2bx],dim=lat).sortby(lat)
			psix2 = psi1x - psi2x
			if method == 'v':
				psi_out = psix2*np.pi/180
			else:
				psi_out = 0.5*(psix1+psix2)*np.pi/180
		# now smooth if required
		if smooth is not None:
			for rdim in smooth.keys():
				psi_out = psi_out.rolling({rdim:smooth[rdim]},center=True,min_periods=1).mean()
		return psi_out

##############################################################################################
def ComputeWaveActivityFlux(phi_or_u,phiref_or_v,uref=None,vref=None,lat='infer',lon='infer',pres='infer',tref=None,is_anomaly=False,qg=False,use_windspharm=False,kwpsi={}):
	'''
		Compute Wave Activity Flux as in Takaya & Nakamura GRL 1997 and Takaya & Nakamura JAS 2001.
		Results checked against plots at http://www.atmos.rcast.u-tokyo.ac.jp/nishii/programs/
		This is 3D wave flux in a non-zonally symmetric base flow, optionally QG approximated.
		Inputs are xarray DataArrays.
		Latitude, longitude in degrees, pressure in hPa.

		INPUTS:
		phi_or_u   : geopotential [m2/s2] (qg = True, typically phi = g*Z) OR zonal wind [m/s] (qg = False)
		phiref_or_v: reference geopotential [m2/s2] (qg = True)
						OR (full) meridional wind [m/s] (qg = False)
		uref	   : reference zonal wind [m/s] - only needed if (qg = False)
		vref	   : reference meridional wind [m/s] - only needed if (qg = False)
		lat	       : name of latitude in DataArrays
		lon	       : name of longitude in DataArrays OR value of pressure level.
		               used for vector scaling on horizontal and for vertical derivatives if 3D
		pres	   : name of pressure in DataArrays [hPa]
		tref	   : reference temperature [K] for static stability parameter S2.
					   If None, only compute horizontal wave flux.
					   Else, used to compute S2.
						 If xr.Dataset or xr.DataArray, compute S2 assuming tref is temperature [K].
						 Else, assume S2 = tref.
						 Note that S2 should be a function of
						  pressure only (see Vallis 2017, Eq 5.127)
                is_anomaly : if True, phi_or_u and phiref_or_v is anomaly. If False, they are full fields. [False]
		qg	       : use QG streamfunction (phi-phiref)/f? Otherwise, use windspharm.streamfunction.
					   Note that if qg=False, phi_or_u = u and phiref_or_v = v to compute the streamfunction.
		use_windspharm: use windspharm package to compute streamfunction. Has no effect if qg=True.
		OUTPUTS:
		Wx,Wy	 : Activity Vectors along lon [m2/s2], lat [m2/s2]
		Wx,Wy,Wz : Activity Vectors along lon [m2/s2], lat [m2/s2], pres [hPa.m/s2]
		div	 : Divergence of horizontal Wave Activity Flux [m/s/day]
		div3	 : Divergence of vertical Wave Activity Flux [m/s/day]
		'''
	import numpy as np
	from xarray import DataArray
	from .constants import a0,Rd,kappa,Omega
	for var in [phi_or_u,uref,vref,phiref_or_v,tref]:
		if not isinstance(var,DataArray):
			if var is None or np.isscalar(var):
				pass
			else:
				# raise ValueError('all inputs have to be xarray.DataArrays!')
				print('ALL INPUTS HAVE TO BE XARRAY.DATAARRYAS! WILL CONTINUE AND TRY MY BEST...')
	# shape
	dim_names = FindCoordNames(phi_or_u)
	if lon == 'infer':
		lon = dim_names['lon']
	if lat == 'infer':
		lat = dim_names['lat']
	if pres == 'infer':
		pres = dim_names['pres']
	p0 = 1.e3
	rad2deg = 180/np.pi
	radlat = np.deg2rad(phi_or_u[lat])
	coslat = np.cos(radlat)
	one_over_coslat2 = coslat**(-2)
	one_over_a2 = a0**(-2)
	one_over_acoslat = (a0*coslat)**(-1)
	f  = 2*Omega*np.sin(radlat)

	if qg:
		use_windspharm = False
	if use_windspharm:
		try:
			from windspharm.xarray import VectorWind
		except:
			use_windspharm = False
			print('WARNING: YOU REQUESTED USE_WINDSPHARM=TRUE, BUT I CANNOT IMPORT WINDSPHARM. CONTINUING WIHOUT WINDSPHARM.')
	if qg:
		uref = -(phiref_or_v/f).differentiate(lat,edge_order=2).reduce(np.nan_to_num)
		uref = uref/a0
		vref =  (phiref_or_v/f).differentiate(lon,edge_order=2).reduce(np.nan_to_num)
		vref = vref*one_over_acoslat
		if is_anomaly:
			psi = phi_or_u/f #[m2/s]
		else:
			psi = (phi_or_u-phiref_or_v)/f #[m2/s]
	else:
		if is_anomaly:
			u = phi_or_u
			v = phiref_or_v
		else:
			u = phi_or_u - uref
			v = phiref_or_v - vref
		psi = ComputeStreamfunction(u,v,lat,lon,use_windspharm=use_windspharm,**kwpsi) #[m2/s]
	# wave activity flux only valid for westerlies.
	#  it is common practice to also mask weak westerlies,
	#  with a value of 1m/s often seen
	mask = uref > 1.0
	mag_u = np.sqrt(uref**2 + vref**2)
	# psi.differentiate(lon) == np.gradient(psi)/np.gradient(lon) [psi/lon]
	# this (with qg=True,use_windspharm=False) has been checked against Fig (a) at
	#   http://www.atmos.rcast.u-tokyo.ac.jp/nishii/programs/
	dpsi_dlon = psi.differentiate(lon,edge_order=2).reduce(np.nan_to_num) # [m2/s/deg]
	dpsi_dlat = psi.differentiate(lat,edge_order=2).reduce(np.nan_to_num) # [m2/s/deg]
	d2psi_dlon2 = dpsi_dlon.differentiate(lon,edge_order=2) # [dpsi_dlon/lon] = [m2/s/deg2]
	d2psi_dlat2 = dpsi_dlat.differentiate(lat,edge_order=2) # [m2/s/deg2]
	d2psi_dlon_dlat = dpsi_dlon.differentiate(lat,edge_order=2) # [m2/s/deg2]

	wx =	uref*one_over_coslat2*one_over_a2*(dpsi_dlon**2 - psi*d2psi_dlon2) \
	    + vref*one_over_a2/coslat*(dpsi_dlon*dpsi_dlat - psi*d2psi_dlon_dlat)
	#[m/s*1/m2*m4/s2/deg2] = [m3/s3/deg2]
	wy =	uref*one_over_a2/coslat*(dpsi_dlon*dpsi_dlat - psi*d2psi_dlon_dlat) \
	    + vref*one_over_a2*(dpsi_dlat**2 - psi*d2psi_dlat2)
	#[m3/s3/deg2]

	# scaling: In TN01, Eq (38) uses spherical coordinates with z in vertical.
	#  the scaling factor is then pref/p0*coslat/2/mag_u.
	#  Eq (C5) shows the same in pressure coordinates, but with (x,y) in the
	#  horizontal. The scaling factor is then 1/2/mag_u.
	# By testing, I found that the latter is consistent with the horizontal EP flux,
	#  i.e. wy.mean('lon') == ep1 if uref = u.mean('lon') if scaling factor = 1/2/mag_u.
	#  Therefore, I use this as scaling factor. Note that this will then cause a
	#   difference compared to Fig (a) at
	#   http://www.atmos.rcast.u-tokyo.ac.jp/nishii/programs/

	coeff = rad2deg**2/2/mag_u # coefficient for p-coords, [deg2*s/m]

	# get the vectors in physical units of m2/s2, correcting for radians vs. degrees
	wx = coeff*wx #[deg2*s/m*m3/s3/deg2] = [m2/s2]
	wy = coeff*wy

	wx.name = 'wx'
	wx.attrs['units'] = 'm2/s2'
	wy.name = 'wy'
	wy.attrs['units'] = 'm2/s2'
	#
	# Now compute the divergence
	#  there is a good chance wx,wy contain NaNs (no propagation where u<0).
	#  therefore, using windspharm for the gradients does not work and we fall
	#  back to the conventional way of computing gradients.
	div1 = wx.differentiate(lon,edge_order=2)
	div2 = (wy*coslat).differentiate(lat,edge_order=2)
	div1 = one_over_acoslat*div1*rad2deg
	div2 = div2*rad2deg/a0

	div = (div1+div2)*86400
	div.name = 'div'
	div.attrs['units'] = 'm/s/d'

	if tref is None:
		return wx.where(mask),wy.where(mask),div.where(mask)
	else:
		# psi.differentiate(pres) == np.gradient(psi)/np.gradient(pres) [psi/pres]
		dpsi_dpres = psi.differentiate(pres,edge_order=2).reduce(np.nan_to_num) # [m2/s/hPa]
		d2psi_dlon_dpres = dpsi_dlon.differentiate(pres,edge_order=2) # [m2/s/deg/hPa]
		d2psi_dlat_dpres = dpsi_dlat.differentiate(pres,edge_order=2) # [m/s/deg/hPa]
		# S2 = -\alpha*\partial_p\ln\theta, \alpha = 1/\rho = Rd*T/p
		#    = R/p*(p/p0)**\kappa d\theta/dp, Vallis (2017, p. 192 (eq. 5.127))
		#  this should be a reference profile and a function of pressure only!
		pressure = uref[pres]
		if isinstance(tref,DataArray):
			pp0 = (p0/pressure)**kappa # []
			theta = pp0*tref # [K]
			S2 = -Rd/pressure*(pressure/p0)**kappa*theta.differentiate(pres,edge_order=2)
			# [J/kg/K/hPa*K/hPa] = [J/kg/hPa2] = [kg*m2/s2/kg/hPa2] = [m2/s2/hPa2] = [m3/kg/hPa]
			# S2 = S2.where(S2>1e-7,1e-7) # avoid division by zero
		else:
			S2 = tref

		wz = f**2/S2*( uref*one_over_acoslat*(dpsi_dlon*dpsi_dpres - psi*d2psi_dlon_dpres) \
					+ vref/a0*(dpsi_dlat*dpsi_dpres - psi*d2psi_dlat_dpres) )
		# units: using [S2] = m3/kg/hPa
		# s-2 * kg * hPa * m-3 * m * s-1 * m-1 * m2 * s-1 * m2 * s-1 * deg-1 * hPa-1 = kg*m/s5/deg

		coeff = rad2deg/2/mag_u # coefficient for p-coords, [deg*s/m]

		wz = coeff*wz # [wz] = kg*m/s5/deg*deg*s/m = kg/s4 = hPa.m/s2.

		wz.name = 'wz'
		wz.attrs['units'] = 'hPa.m/s2'
		wz.attrs['alternative_units'] = 'kg/s4'

		div3 = wz.differentiate(pres,edge_order=2)*86400
		div3.name = 'div'
		div3.attrs['units'] = 'm/s/d'

		return wx.where(mask),wy.where(mask),wz.where(mask),div.where(mask),div3.where(mask)

##############################################################################################
def PlotEPfluxArrows(x,y,ep1,ep2,fig,ax,xlim=None,ylim=None,xscale='linear',yscale='linear',invert_y=True, newax=False, pivot='tail',scale=None,quiv_args=None):
	"""Correctly scales the Eliassen-Palm flux vectors for plotting on a latitude-pressure or latitude-height axis.
		x,y,ep1,ep2 assumed to be xarray.DataArrays.

	INPUTS:
		x	: horizontal coordinate, assumed in degrees (latitude) [degrees]
		y	: vertical coordinate, any units, but usually this is pressure or height
		ep1	: horizontal Eliassen-Palm flux component, in [m2/s2]. Typically, this is ep1_cart from
				   ComputeEPfluxDiv()
		ep2	: vertical Eliassen-Palm flux component, in [U.m/s2], where U is the unit of y.
				   Typically, this is ep2_cart from ComputeEPfluxDiv(), in [hPa.m/s2] and y is pressure [hPa].
		fig	: a matplotlib figure object. This figure contains the axes ax.
		ax	: a matplotlib axes object. This is where the arrows will be plotted onto.
		xlim	: axes limits in x-direction. If None, use [min(x),max(x)]. [None]
		ylim	: axes limits in y-direction. If None, use [min(y),max(y)]. [None]
		xscale	: x-axis scaling. currently only 'linear' is supported. ['linear']
		yscale	: y-axis scaling. 'linear' or 'log' ['linear']
		invert_y: invert y-axis (for pressure coordinates). [True]
		newax	: plot on second y-axis. [False]
		pivot	: keyword argument for quiver() ['tail']
		scale	: keyword argument for quiver(). Smaller is longer [None]
				  besides fixing the length, it is also usefull when calling this function inside a
				   script without display as the only way to have a quiverkey on the plot.
               quiv_args: further arguments passed to quiver plot.

	OUTPUTS:
	   Fphi*dx : x-component of properly scaled arrows. Units of [m3.inches]
	   Fp*dy   : y-component of properly scaled arrows. Units of [m3.inches]
	   ax	: secondary y-axis if newax == True
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	#
	def Deltas(z,zlim):
		# if zlim is None:
		return np.max(z)-np.min(z)
		# else:
			# return zlim[1]-zlim[0]
	# Scale EP vector components as in Edmon, Hoskins & McIntyre JAS 1980:
	cosphi = np.cos(np.deg2rad(x))
	a0 = 6376000.0 # Earth radius [m]
	grav = 9.81
	# first scaling: Edmon et al (1980), Eqs. 3.1 & 3.13
	Fphi = 2*np.pi/grav*cosphi**2*a0**2*ep1 # [m3.rad]
	Fp   = 2*np.pi/grav*cosphi**2*a0**3*ep2 # [m3.hPa]
	#
	# Now comes what Edmon et al call "distances occupied by 1 radian of
	#  latitude and 1 [hecto]pascal of pressure on the diagram."
	# These distances depend on figure aspect ratio and axis scale
	#
	# first, get the axis width and height for
	#  correct aspect ratio
	width,height = GetAxSize(fig,ax)
	# we use min(),max(), but note that if the actual axis limits
	#  are different, this will be slightly wrong.
	delta_x = Deltas(x,xlim)
	delta_y = Deltas(y,ylim)
	#
	#scale the x-axis:
	if xscale == 'linear':
		dx = width/delta_x/np.pi*180
	else:
		raise ValueError('ONLY LINEAR X-AXIS IS SUPPORTED AT THE MOMENT')
	#scale the y-axis:
	if invert_y:
		y_sign = -1
	else:
		y_sign = 1
	if yscale == 'linear':
		dy = y_sign*height/delta_y
	elif yscale == 'log':
		dy = y_sign*height/y/np.log(np.max(y)/np.min(y))
	#
	# plot the arrows onto axis
	quivArgs = {'angles':'uv','scale_units':'inches','pivot':pivot}
	if quiv_args is not None:
		for key in quiv_args.keys():
			quivArgs[key] = quiv_args[key]
	if scale is not None:
		quivArgs['scale'] = scale
	if newax:
		ax = ax.twinx()
		ax.set_ylabel('pressure [hPa]')
	try:
		Q = ax.quiver(x,y,Fphi*dx,Fp*dy,**quivArgs)
	except:
		Q = ax.quiver(x,y,dx*Fphi.transpose(),dy*Fp.transpose(),**quivArgs)
	if scale is None:
		fig.canvas.draw() # need to update the plot to get the Q.scale
		U = Q.scale
	else:
		U = scale
	if U is not None: # when running inside a script, the figure might not exist and therefore U is None
		ax.quiverkey(Q,0.9,1.02,U/width,label=r'{0:.1e}$\,m^3$'.format(U),labelpos='E',coordinates='axes')
	if invert_y:
		ax.invert_yaxis()
	if xlim is not None:
		ax.set_xlim(xlim)
	if ylim is not None:
		ax.set_ylim(ylim)
	ax.set_yscale(yscale)
	ax.set_xscale(xscale)
	#
	if newax:
		return Fphi*dx,Fp*dy,ax
	else:
		return Fphi*dx,Fp*dy

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
	  integ- averaged data, length N
	"""
	#make sure there are more than one grid points
	if len(lat) < 2:
		return np.mean(data,axis=axis)
	#get data into the correct shape
	tmp = AxRoll(data,axis)
	shpe= tmp.shape
	tmp = np.reshape(tmp,(shpe[0],np.prod(shpe[1:])))
	#cosine weighting
	J = np.where((lat>=lim)*(lat<=mx))[0]
	coslat = np.cos(np.deg2rad(lat))**cosp
	coswgt = np.trapz(coslat[J],lat[J])
	tmp = np.trapz(tmp[J,:]*coslat[J][:,np.newaxis],lat[J],axis=0)/coswgt
	integ = np.reshape(tmp,shpe[1:])
	return integ

##############################################################################################
def GlobalAvgXr(ds,lats=[-90,90],lat='infer',cosp=1):
	"""Compute cosine weighted meridional average.

	INPUTS:
	  ds	- xarray.Dataset or xarray.DataArray to average.
	  lats	- latitude range to average over.
	  lat	- name of latitude dimension - guess if 'infer'.
	  cosp	- power of cosine weighting
	OUTPUTS:
	  averaged data, size ds minus latitude
	"""
	# try to find latitude dimension
	if lat == 'infer':
		lat = FindCoordNames(ds)['lat']
	# make sure order of latitude is correct
	revert = False
	if ds[lat][0] < ds[lat][-1]:
		if lats[0] > lats[1]:
			revert = True
	else:
		if lats[0] < lats[1]:
			revert = True
	if revert:
		lats = [lats[1],lats[0]]
	dt = ds.sel({lat:slice(*lats)})
	coslat = np.cos(np.deg2rad(dt[lat]))**cosp
	dt.attrs['average'] = 'cos^{2} weighted meridional average from {0} to {1}'.format(lats[0],lats[1],cosp)
	return dt.weighted(coslat).mean(lat,keep_attrs=True)

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
	dp   = np.gradient(pres)[:,np.newaxis]*100.
	dTdp = np.gradient(Tz,edge_order=2)[0]/dp
	p = pres[:,np.newaxis]*100. # [Pa]
	N2 = -Rd*p/(H**2.) * (dTdp - Rd*Tz/(p*cp))
	return N2

##############################################################################################
def ComputeN2Xr(Tz,pres='infer',H=7.e3,Rd=287.04,cp=1004):
	''' Compute the Brunt-Vaisala frequency from zonal mean temperature
		 N2 = -Rd*p/(H**2.) * (dTdp - Rd*Tz/(p*cp))
		 this is equivalent to
		 N2 = g/\theta d\theta/dz, with p = p0 exp(-z/H)

		INPUTS:
			Tz    - zonal mean temperature [K], xarray.DataArray
			pres  - name of pressure [hPa]
			H     - scale height [m]
			Rd    - specific gas constant for dry air
			cp    - specific heat of air at constant pressure
		OUTPUTS:
			N2  - Brunt-Vaisala frequency, [1/s2], dim pres x lat
	'''
	if pres == 'infer':
		dim_names = FindCoordNames(Tz)
		pres = dim_names['pres']
	dTdp = Tz.differentiate(pres,edge_order=2)*0.01
	p = Tz[pres]*100. # [Pa]
	N2 = -Rd*p/(H**2.) * (dTdp - Rd*Tz/(p*cp))
	return N2

##############################################################################################
# @jit
def FlexiGradPhi(data,dphi):
	if len(data.shape) == 3:
		grad = np.gradient(data,edge_order=2)[2]
	else:
		grad = np.gradient(data,edge_order=2)[1]
	return grad/dphi
# @jit
def FlexiGradP(data,dp):
	if len(data.shape) == 3:
		grad = np.gradient(data,edge_order=2)[1]
	else:
		grad = np.gradient(data,edge_order=2)[0]
	return grad/dp
# @jit(cache=True)
def ComputeMeridionalPVGrad(lat, pres, uz, Tz, Rd=287.04, cp=1004, a0=6.371e6, component='ABC'):
	'''Compute the meridional gradient of potential vorticity.
		Computed following Simpson et al JAS (2009) DOI 10.1175/2008JAS2758.1.
		This quantity has three terms,
		q_\phi = A - B + C, where
				A = 2*Omega*cos\phi
				B = \partial_\phi[\partial_\phi(ucos\phi)/acos\phi]
				C = af^2/Rd*\partial_p(p\theta\partial_pu/(T\partial_p\theta))

		INPUTS:
			lat	  - latitude [degrees]
			pres	  - pressure [hPa]
			uz	  - zonal mean zonal wind [m/s], dim pres x lat OR N x pres x lat
			Tz	  - zonal mean temperature [K], dim pres x lat OR N x pres x lat
			component - option to only return one, two, or all of the components.
						 Add a letter for each of the components 'A', 'B', 'C'.
						 Note: As B has a minus sign in q_\phi, option 'B' returns -B
		OUTPUTS:
			q_phi - meridional gradient of potential vorticity [1/s], dim pres x lat OR N x pres x lat
	'''
	if not ('A' in component)+('B' in component)+('C' in component):
		raise ValueError('component has to contain A,B and/or C, but got '+component)
	# some constants
	from .constants import Omega
	from .constants import p0_Pa as p0

	## make sure we have the dimesions as expected
	if uz.shape != Tz.shape:
		raise ValueError('UZ AND TZ DO NOT HAVE THE SAME SHAPE')
	elif len(uz.shape) > 3:
		raise ValueError('TOO MANY DIMENSIONS IN UZ AND TZ')

	## convert to Pa
	p = pres[:]*100
	if len(uz.shape) == 3:
		dp = np.gradient(p)[np.newaxis,:,np.newaxis]
		p  = p[np.newaxis,:,np.newaxis]
	else:
		dp = np.gradient(p)[:,np.newaxis]
		p = p[:,np.newaxis]
	## convert to radians
	latpi = np.deg2rad(lat)
	if len(uz.shape) == 3:
		dphi = np.gradient(latpi)[np.newaxis,np.newaxis,:]
		latpi= latpi[np.newaxis,np.newaxis,:]
	else:
		dphi = np.gradient(latpi)[np.newaxis,:]
		latpi = latpi[np.newaxis,:]

	#
	result = np.zeros(uz.shape)
	## first term A
	if 'A' in component:
		A = 2*Omega*np.cos(latpi)
		result += A

	#
	## second term B
	if 'B' in component:
		dudphi = FlexiGradPhi(uz*np.cos(latpi),dphi)
		B = dudphi/np.cos(latpi)/a0
		B = FlexiGradPhi(B,dphi)
		result -= B

	#
	## third term C
	if 'C' in component:
		f = 2*Omega*np.sin(latpi)

		dudp = FlexiGradP(uz,dp)

		kappa = Rd/cp
		pp0   = (p0/p)**kappa
		theta = Tz*pp0
		theta_p = FlexiGradP(theta,dp)

		C = p*theta*dudp/(Tz*theta_p)
		C = FlexiGradP(C,dp)
		C = a0*f*f*C/Rd
		result += C
	#
	return result

##############################################################################################
def ComputeMeridionalPVGradXr(uz, Tz, lat='lat', pres='level', Rd=287.04, cp=1004, a0=6.371e6, component='ABC'):
	'''Compute the meridional gradient of potential vorticity.
		Computed following Simpson et al JAS (2009) DOI 10.1175/2008JAS2758.1.
		This quantity has three terms,
		q_\phi = A - B + C, where
				A = 2*Omega*cos\phi
				B = \partial_\phi[\partial_\phi(ucos\phi)/acos\phi]
				C = af^2/Rd*\partial_p(p\theta\partial_pu/(T\partial_p\theta))

		INPUTS:
			uz        - zonal mean zonal wind [m/s], dim pres x lat OR N x pres x lat
			Tz        - zonal mean temperature [K], dim pres x lat OR N x pres x lat
			lat       - name of latitude [degrees]
			pres      - name of pressure [hPa]
			component - option to only return one, two, or all of the components.
						 Add a letter for each of the components 'A', 'B', 'C'.
						 Note: As B has a minus sign in q_\phi, option 'B' returns -B
		OUTPUTS:
			q_phi - meridional gradient of potential vorticity [1/s], dim pres x lat OR N x pres x lat
	'''
	if not ('A' in component)+('B' in component)+('C' in component):
		raise ValueError('component has to contain A,B and/or C, but got '+component)
	# some constants
	from .constants import Omega,p0
	coslat = np.cos(np.deg2rad(uz[lat]))
	factor_phi = 180/np.pi
	factor_pres= 0.01

	#
	result = uz*0
	## first term A
	if 'A' in component:
		A = 2*Omega*coslat
		result += A

	#
	## second term B
	if 'B' in component:
		dudphi = (coslat*uz).differentiate(lat,edge_order=2)*factor_phi
		B = dudphi/coslat/a0
		B = B.differentiate(lat,edge_order=2)*factor_phi
		result -= B

	#
	## third term C
	if 'C' in component:
		f = 2*Omega*np.sin(np.deg2rad(uz[lat]))

		dudp = uz.differentiate(pres,edge_order=2)*factor_pres

		kappa = Rd/cp
		pp0   = (p0/uz[pres])**kappa
		theta = Tz*pp0
		theta_p = theta.differentiate(pres,edge_order=2)*factor_pres

		C = uz[pres]*theta*dudp/(Tz*theta_p)
		C = C.differentiate(pres,edge_order=2)*factor_pres
		C = a0*f*f*C/Rd
		result += C
	#
	return result


##############################################################################################
def ComputeRefractiveIndex(lat,pres,uz,Tz,k,N2const=None):
	'''
		Refractive index as in Simpson et al (2009) doi 10.1175/2008JAS2758.1 and also Matsuno (1970) doi 10.1175/1520-0469(1970)027<0871:VPOSPW>2.0.CO;2
		Stationary waves are assumed, ie c=0.

		Setting k=0 means the only term depending on wave number is left out. This could be more efficient if n2(k) for different values of k is of interest.

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
			k     - zonal wave number. [.]
			N2const - if not None, assume N2 = const = N2const [1/s2]
		Outputs are:
			n2  - refractive index, dimension pres x lat [.]
	'''
	# some constants
	from .constants import Rd,cp,a0,Omega
	H     = 7.e3 # [m]

	latpi = np.deg2rad(lat)

	#
	## term D
	dqdy = ComputeMeridionalPVGrad(lat,pres,uz,Tz,Rd,cp,a0)
	D = dqdy/(a0*uz)

	#
	## term E
	latpi = latpi[np.newaxis,:]
	E = ( k/(a0*np.cos(latpi)) )**2

	#
	## term F
	f = 2*Omega*np.sin(latpi)
	f2 = f*f
	if N2const is None:
		N2 = ComputeN2(pres,Tz,H,Rd,cp)
	else:
		N2 = N2const
	H2 = H*H
	F = f2/(4*N2*H2)

	return a0*a0*(D-E-F)

##############################################################################################
def ComputeRefractiveIndexXr(uz,Tz,k,lat='lat',pres='level',N2const=None,ulim=None):
	'''
		Refractive index as in Simpson et al (2009) doi 10.1175/2008JAS2758.1 and also Matsuno (1970) doi 10.1175/1520-0469(1970)027<0871:VPOSPW>2.0.CO;2
		Stationary waves are assumed, ie c=0.

		Setting k=0 means the only term depending on wave number is left out. This could be more efficient if n2(k) for different values of k is of interest.

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
			uz    - zonal mean zonal wind, xarray.DataArray
			Tz    - zonal mean temperature, xarray.DataArray
			k     - zonal wave number. [.]
			lat   - name of latitude [degrees]
			pres  - name of pressure [hPa]
			N2const - if not None, assume N2 = const = N2const [1/s2]
			ulim  - only compute n2 where |u| > ulim to avoid divisions by zero.
		Outputs are:
			n2  - refractive index, dimension pres x lat [.]
	'''
	# some constants
	from .constants import Rd,cp,a0,Omega
	H     = 7.e3 # [m]

	#
	## term D
	dqdy = ComputeMeridionalPVGradXr(uz,Tz,lat,pres,Rd,cp,a0)
	if ulim is not None:
		utmp = uz.where(np.abs(uz)>ulim)
	else:
		utmp = uz
	D = dqdy/(a0*utmp)

	#
	## term E
	coslat = np.cos(np.deg2rad(uz[lat]))
	E = ( k/(a0*coslat) )**2

	#
	## term F
	sinlat = np.sin(np.deg2rad(uz[lat]))
	f = 2*Omega*sinlat
	f2 = f*f
	if N2const is None:
		N2 = ComputeN2Xr(Tz,pres,H,Rd,cp)
	else:
		N2 = N2const
	H2 = H*H
	F = f2/(4*N2*H2)

	return a0*a0*(D-E-F)

##############################################################################################
def GetWaves(x,y=None,wave=-1,axis=-1,do_anomaly=False):
	"""Get Fourier mode decomposition of x, or <x*y>, where <.> is zonal mean.

		If y!=[], returns Fourier mode contributions (amplitudes) to co-spectrum zonal mean of x*y. Shape is same as input, except axis which is len(axis)/2+1 due to Fourier symmetry for real signals.

		If y=[] and wave>=0, returns real space contribution of given wave mode. Output has same shape as input.
		If y=[] and wave=-1, returns real space contributions for all waves. Output has additional first dimension corresponding to each wave.

	INPUTS:
		x	   - the array to decompose
		y	   - second array if wanted
		wave	   - which mode to extract. all if <0
		axis	   - along which axis of x (and y) to decompose
		do_anomaly - decompose from anomalies or full data
	OUTPUTS:
		xym	   - data in Fourier space
	"""
	initShape = x.shape
	x = AxRoll(x,axis)
	if y is not None:
		y = AxRoll(y,axis)
	# compute anomalies
	if do_anomaly:
		x = GetAnomaly(x,0)
		if y is not None:
			y = GetAnomaly(y,0)
	# Fourier decompose
	x = np.fft.fft(x,axis=0)
	nmodes = x.shape[0]//2+1
	if wave < 0:
			if y is not None:
				xym = np.zeros((nmodes,)+x.shape[1:])
			else:
				xym = np.zeros((nmodes,)+initShape)
	else:
		xym = np.zeros(initShape[:-1])
	if y is not None:
			y = np.fft.fft(y,axis=0)
			# Take out the waves
			nl  = x.shape[0]**2
			xyf  = np.real(x*y.conj())/nl
			# due to symmetric spectrum, there's a factor of 2, but not for wave-0
			mask = np.zeros_like(xyf)
			if wave < 0:
				for m in range(xym.shape[0]):
					mask[m,:] = 1
					mask[-m,:]= 1
					xym[m,:] = np.sum(xyf*mask,axis=0)
					mask[:] = 0
				# wavenumber 0 is total of all waves
                                #  this makes more sense than the product of the zonal means
				xym[0,:] = np.nansum(xym[1:,:],axis=0)
				xym = AxRoll(xym,axis,invert=True)
			else:
				xym = xyf[wave,:]
				if wave >= 0:
					xym = xym + xyf[-wave,:]
	else:
			mask = np.zeros_like(x)
			if wave >= 0:
				mask[wave,:] = 1
				mask[-wave,:]= 1 # symmetric spectrum for real signals
				xym = np.real(np.fft.ifft(x*mask,axis=0))
				xym = AxRoll(xym,axis,invert=True)
			else:
				for m in range(xym.shape[0]):
					mask[m,:] = 1
					mask[-m,:]= 1 # symmetric spectrum for real signals
					fourTmp = np.real(np.fft.ifft(x*mask,axis=0))
					xym[m,:] = AxRoll(fourTmp,axis,invert=True)
					mask[:] = 0
	return np.squeeze(xym)

##helper functions
def GetAnomaly(x,axis=-1):
	"""Computes the anomaly of array x along dimension axis.

	INPUTS:
	  x    - array to compute anomalies from
	  axis - axis along dimension for anomalies
	OUTPUTS:
	  x    - anomalous array
	"""    #bring axis to the front
	xt= AxRoll(x,axis)
	#compute anomalies
	xt = xt - xt.mean(axis=0)[np.newaxis,:]
	#bring axis back to where it was
	x = AxRoll(xt,axis,invert=True)
	return x

##############################################################################################
def GetWavesXr(x,y=None,wave=-1,dim='infer',anomaly=None):
	"""Get Fourier mode decomposition of x, or <x*y>, where <.> is zonal mean.

		If y!=None, returns Fourier mode contributions (amplitudes) to co-spectrum zonal mean of x*y. Dimension along which Fourier is performed is either gone (wave>=0) or has len(axis)/2+1 due to Fourier symmetry for real signals (wave<0).

		If y=None and wave>=0, returns real space contribution of given wave mode. Output has same shape as input.
		If y=None and wave<0, returns real space contributions for all waves. Output has additional first dimension corresponding to each wave.

	INPUTS:
		x	   - the array to decompose. xr.DataArray
		y	   - second array if wanted. xr.DataArray
		wave	   - which mode to extract. all if <0
		dim	   - along which dimension of x (and y) to decompose. Defaults to longitude if 'infer'.
		anomaly	   - if not None, name of dimension along which to compute anomaly first.
	OUTPUTS:
		xym	   - data. xr.DataArray
	"""
	from xarray import DataArray
	if dim == 'infer':
		dim_names = FindCoordNames(x)
		dim = dim_names['lon']
	if anomaly is not None:
		x = x - x.mean(anomaly)
		if y is not None:
			y = y - y.mean(anomaly)
	sdims = [d for d in x.dims if d != dim]
	if len(sdims) == 0:
		xstack = x.expand_dims('stacked',axis=-1)
	else:
		xstack = x.stack(stacked=sdims)
	if y is None:
		ystack=None
	else:
		if len(sdims) == 0:
			ystack = y.expand_dims('stacked',axis=-1)
		else:
			ystack = y.stack(stacked=sdims)
	gw = GetWaves(xstack,ystack,wave=wave,axis=xstack.get_axis_num(dim))
	if y is None and wave >= 0: # result in real space
		if len(sdims) == 0:
			stackcoords = x.coords
		else:
			stackcoords = [xstack[d] for d in xstack.dims]
	elif y is None and wave < 0: # additional first dimension of wave number
		stackcoords = [('k',np.arange(gw.shape[0]))]
		if len(sdims) == 0:
			stackcoords = stackcoords + [x[d] for d in x.dims]
		else:
			stackcoords = stackcoords + [xstack[d] for d in xstack.dims]
	elif y is not None and wave >= 0: # original dimension is gone
		stackcoords = [xstack.stacked]
	elif y is not None and wave < 0: # additional dimension of wavenumber
		stackcoords = [('k',np.arange(gw.shape[0])), xstack.stacked]
	gwx = DataArray(gw,coords=stackcoords)
	return gwx.unstack()



#######################################################
def Meters2Coord(data,coord,mode='m2lat',axis=-1):
	"""Convert value (probably vector component) of one unit
		into another one.

		INPUTS:
		data  - data to convert
		coord - values of latitude [degrees] or pressure [hPa]
		mode  - 'm2lat', 'm2lon', 'm2hPa', and all inverses
		axis  - axis of data which needs modification
		OUTPUTS:
		out   - converted from data
		"""
	# constants
	from .constants import a0,p0
	ps    = p0
	H     = 7e3
	# geometric quantitites
	if 'lat' in mode or 'lon' in mode:
	   rad2deg = 180/np.pi
	   coslat  = np.cos(coord/rad2deg)
	   cosm1   = 1/coslat
	   gemfac  = rad2deg/a0
	#
	ndims = len(np.shape(data))
	if mode == 'm2lat':
		out = data*gemfac
	elif mode == 'lat2m':
		out = data/gemfac
	elif mode in ['m2lon','lon2m','m2hPa','hPa2m']:
		if ndims > 1:
			tmp = AxRoll(data,axis)
	# else:
	#	tmp = data
	#	out = np.zeros_like(tmp)
	else:
		raise ValueError("mode not recognized")
	if mode == 'm2lon':
		if ndims > 1:
			for l in range(out.shape[0]):
				out[l,:] = tmp[l,:]*cosm1
		else:
			out = tmp*cosm1
		out = out*gemfac
		out = AxRoll(out,axis,invert=True)
	elif mode == 'lon2m':
		if ndims > 1:
			for l in range(out.shape[0]):
				out[l,:] = tmp[l,:]*coslat
		else:
			out = tmp*coslat
		out = out/gemfac
	elif mode == 'm2hPa':
		if ndims > 1:
			for p in range(out.shape[0]):
				out[p,:] = -coord[p]*tmp[p,:]
		else:
			out = -coord*tmp
		out = out/H
		out = AxRoll(out,axis,invert=True)
	elif mode == 'hPa2m':
		if ndims > 1:
			for p in range(out.shape[0]):
				out[p,:] = -coord[p]/tmp[p,:]
		else:
			out = -coord/tmp
		out[tmp==0] = NaN
		out = out*H
		out = AxRoll(out,axis,invert=True)
	#
	return out

##############################################################################################
def ComputeBaroclinicity(lat, tempIn, hemi='both', minLat=20, maxLat=60, pres=None, minPres=250):
	"""Compute the meridional temperature gradient, integrated between minLat and maxLat.
		Thus, baroclinicity is defined as the difference in temperature between minLat and maxLat,
		optionally integrated from the surface to minPres.

		INPUTS:
		lat	- latitude [degrees]
		tempIn	- temperature, shape time x pressure x lat
		hemi	- hemisphere to consider. 'both','N','S'
		minLat	- latitude closest to equator
		maxLat	- latitdue closest to pole
		pres	- pressure [hPa]
		minPres - top of pressure averaging from surface.
					only used if pres is not None.
		OUTPUT:
		dT	- dictionary with options ['S'] and ['N'] if hemi='both'
					otherwise array of length len(time)
	"""

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
		Convert an array of times with given units and calendar
		 into the same times but in units of days.
	"""
	import netCDF4 as nc
	import netcdftime as nct
	date = nc.num2date(time,units,calendar)
	unitArray = units.split()
	dayUnits = units.replace(unitArray[0],'days')
	t = nct.utime(dayUnits,calendar=calendar)
	return t.date2num(date)

#######################################################
def StandardGrid(data,lon_name='infer',lat_name='infer',pres_name='infer',rename=False,pdir=None):
	"""
		Make sure longitude is in [0,360] and latitude sorted
		 from lowest to highest.
		 Assumes an xr.DataArray or xr.Dataset.

		INPUTS:
			data:	   xarray.DataArray or xarray.Dataset to regrid
			lon_name:  name of longitude dimension.
						If None: nothing should be done.
						If 'infer': try to find name of longitude
			lat_name:  name of latitude dimension.
						If None: nothing should be done.
						If 'infer': try to find name of latitude
			pres_name: name of pressure dimension.
						If None: nothing should be done.
						If 'infer': try to find name of pressure
			rename:	   Rename dims to standard? If True, output names will be:
						longitude = 'lon'
						latitude  = 'lat'
						pressure  = 'press'
			pdir:	   Re-arrange pressure (or other vertical coordinate):
						If None: don't do anything
						If 'decrease': first value is largest
						If 'increase': first value is smallest
		OUTPUTS:
			data:	  xarray.DataArray with latitude from lowest to highest and
					   longitude between 0 and 360 degrees.
	"""
	if lon_name == 'infer' or lat_name == 'infer' or pres_name == 'infer':
		dim_names = FindCoordNames(data)
		if lon_name == 'infer' and 'lon' in dim_names.keys():
			lon_name = dim_names['lon']
		if lat_name == 'infer' and 'lat' in dim_names.keys():
			lat_name = dim_names['lat']
		if pres_name == 'infer' and 'pres' in dim_names.keys():
			pres_name = dim_names['pres']
	if lat_name is not None and lat_name in data.coords:
		if data[lat_name][0] > data[lat_name][-1]:
			data = data.sortby(lat_name)
	if lon_name is not None and lon_name in data.coords and data[lon_name].min() < 0:
		data = data.assign_coords({lon_name : (data[lon_name]+360)%360})
		data = data.sortby(lon_name)
	if pdir is not None:
		p_name = dim_names['pres']
		if pdir == 'decrease':
			data = data.sortby(p_name,ascending=False)
		elif pdir == 'increase':
			data = data.sortby(p_name,ascending=True)
	if rename:
		data = data.rename(dict(map(reversed, dim_names.items())))
	return data

#######################################################
def ERA2Model(data,lon_name='longitude',lat_name='latitude'):
	"""This function is deprecated. Please see StandardGrid().
	"""
	import warnings
	warnings.warn('ERA2Model() is deprecated. Please use StandardGrid() instead.')
	return StandardGrid(data,lon_name,lat_name)

#######################################################
def ComputeGeostrophicWind(Z,lon_name='longitude',lat_name='latitude',qg_limit=5.0):
	"""
		Compute u,v from geostrophic balance
		Z is geopotential height
	"""
	import numpy as np
	a0 = 6376.0e3
	Omega = 7.292e-5
	grav = 9.81
	# f = 0 at the equator, which will obviously be a problem
	#  so we make f=infty there so that u,v = 0 between -qg_limit and +qg_limit
	f = 2*Omega*np.sin(np.deg2rad(Z[lat_name])).where(np.abs(Z[lat_name])>qg_limit,np.infty)
	cosPhi = np.cos(np.deg2rad(Z[lat_name]))
	u = -grav*Z.differentiate(lat_name)/f/a0
	v =  grav*Z.differentiate(lon_name)/f/a0/cosPhi
	return u,v

#######################################################
def ComputeRossbyWaveSource(u,v):
	"""
		Compute the Rossby Wave Source.
		This is directly from the windspharm documentation,
		https://ajdawson.github.io/windspharm/latest/examples/rws_xarray.html

		INPUTS:
			u : zonal wind [m/s]
			v : meridional wind [m/s]

		OUTPUTS:
			S : Rossby Wave Source [1/s2]
	"""
	from windspharm.xarray import VectorWind
	# Create a VectorWind instance to handle the computations.
	w = VectorWind(u, v)

	# Compute components of rossby wave source: absolute vorticity, divergence,
	# irrotational (divergent) wind components, gradients of absolute vorticity.
	# divergence can contain spurious features ("stripes") at high res. truncate to avoid this.
	dim_names = FindCoordNames(u)
	npnts = min(len(u[dim_names['lat']]),len(u[dim_names['lon']]))
	eta = w.absolutevorticity()
	div = w.divergence(truncation=npnts//2)
	uchi, vchi = w.irrotationalcomponent()
	etax, etay = w.gradient(eta)

	# Combine the components to form the Rossby wave source term.
	rws = eta * -1. * div - (uchi * etax + vchi * etay)
	rws.name = 'rws'
	rws.attrs['units'] = '1/s2'
	rws.attrs['long_name'] = 'Rossby Wave Source'
	return rws

#######################################################
def Projection(projection='EqualEarth',nrows=1,ncols=1,transform='PlateCarree',coast=False,kw_args=None,fig_args=None):
	"""
		Create a new figure with a given projection. To plot data, invoke:
			fig,ax,kw = Projection()
			ax.contourf(x,y,z,**kw)

		INPUTS:
			projection:  desired map projection for plotting.
			nrows:	     how many rows in pyplot.subplots()
			ncols:	     how many columns in pyplot.subplots()
			transform:   map projection of the data.
			coast:	     whether or not to draw coastlines (can be done later with ax.coastines())
			kw_args:     keyword arguments (dict()) for projection.

		OUTPUTS:
			fig:	     figure object
			ax:	     axis with map projection
			transf:	     transform in form of dictionary.
	"""
	from cartopy import crs as ccrs
	cproj = getattr(ccrs,projection)
	if kw_args is None:
		proj = cproj()
	else:
		proj = cproj(**kw_args)
	from matplotlib import pyplot as plt
	if fig_args is None:
		fig,ax = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection':proj})
	else:
		fig,ax = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection':proj},**fig_args)
	if coast:
		if nrows*ncols > 1:
			for a in ax.flatten():
				a.coastlines()
		else:
			ax.coastlines()
	return fig,ax,{'transform':getattr(ccrs,transform)()}

#######################################################
def Cart2Sphere(u, v, w, lon='infer', lat='infer', pres=None, H=7e3, p0=1e3):
	"""
		Convert 3D vector (u,v,w) from Cartesian coordinates to
			spherical coordinates for correct plotting on a 3D sphere.
		Also converts pressure into height if pres is a string with the
			name of the pressure coordinate. Note that this will also transform
			the vertical component w.
		u,v,w must be xarray.DataArrays, and longitude, latitude in degrees.

		INPUTS:
			u,v,w	: x-, y-, and z-components of Cartesian vector
			lon,lat : names of longitude and latitude in DataArray u.
						Angles are in degrees.
			pres	: None: don't do anything with the vertical coodinate
					  string: name of pressure coordinate. Will then be converted
							  to log-pressure with a scale height of H and a base
							  pressure of p0.
			H	: Scale height. A priori, this is in meters. But this parameter
						can be used to adjust for aspect ratio. For instance, if one
						wants more detail in the vertical, making H 10x larger results
						in a 10x "deeper" atmosphere with respect to Earth's radius.
			p0	: Base pressure. This should probably never be changed.

	"""
	from xarray import DataArray
	for var in [u,v,w]:
		if not isinstance(var,DataArray):
			raise ValueError('u,v,w must be xarray DataArrays!')
	import numpy as np
	dim_names = FindCoordNames(u)
	if lon == 'infer':
		lon = dim_names['lon']
	if lat == 'infer':
		lat = dim_names['lat']
	if pres == 'infer':
		pres = dim_names['pres']
	radlat = np.deg2rad(u[lat])
	radlon = np.deg2rad(u[lon])
	coslat = np.cos(radlat)
	sinlat = np.sin(radlat)
	coslon = np.cos(radlon)
	sinlon = np.sin(radlon)
	# pressure to height
	if pres is not None:
		c = H*np.log(u[pres]/(u[pres]+w))
	else:
		c = w
	# transform vector into spherical geometry
	a = -u*sinlon -v*sinlat*coslon +c*coslat*coslon
	b =  u*coslon -v*sinlat*sinlon +c*coslat*sinlon
	c =	       v*coslat	       +c*sinlat
	a.name = u.name
	b.name = v.name
	c.name = w.name
	if pres is None:
		return a,b,c
	# also change the vertical coordinate from pressure to height (log-pressure)
	z = -H*np.log(u[pres]/p0).values
	a = a.assign_coords({pres:z})
	b = b.assign_coords({pres:z})
	c = c.assign_coords({pres:z})
	a = a.rename({pres:'z'})
	b = b.rename({pres:'z'})
	c = c.rename({pres:'z'})
	return a,b,c


#######################################################
def OceanIndex(sst, index='3.4', lon='infer', lat='infer', time='time', avg=5):
	"""
		Produce variability index timeseries from SST according to Technical Notes
		 guidance from UCAR: https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni

		INPUTS:
			sst:  xarray.DataArray which will be averaged over Nino domains
			lon:  name of longitude dimension. Has to be in [0,360]. If 'infer', try to guess.
			lat:  name of latitude dimension. Has to be increasing. If 'infer', try to guess.
			time: name of time dimension. If 'infer', try to guess.
			avg:  size of rolling window for rolling time average [months].
			index: which index to compute. Choices are
				   Nino: '1+2','3','4','3.4','oni','tni','modiki'
					'oni' is the same as '3.4', as they only
					 distinguish themselves via the rolling
					 average period (3 months for oni,
					 5 months for 3.4)
                                        'modiki' is the same as 'tni'
                                    IOD: 'dmi' is Indian Ocean Dipole

		OUTPUTS:
			sst: spatially averaged over respective index domain
				  note that output is monthly
	"""
	# shape
	if lon == 'infer' or lat == 'infer':
		dim_names = FindCoordNames(sst)
	if lon == 'infer':
		lon = dim_names['lon']
	if lat == 'infer':
		lat = dim_names['lat']
	indList = {
		'1+2' : {lon:slice(270,280),lat:slice(-10,0)},
		'3'   : {lon:slice(210,270),lat:slice(-5,5)},
		'4'   : {lon:slice(160,210),lat:slice(-5,5)},
		'3.4' : {lon:slice(190,240),lat:slice(-5,5)},
		'oni' : {lon:slice(190,240),lat:slice(-5,5)},
                'tni' : '1+2std - 4std',
                'modiki':'tni',
                'dmi' : 'dmi1 - dmi2',
                'dmi1': {lon:slice(50,70)  ,lat:slice(-10,10)},
                'dmi2': {lon:slice(90,110) ,lat:slice(-10,0)}
	}
	possible_ninos = list(indList.keys())
	if index not in possible_ninos:
		raise ValueError('Nino type {0} not recognised. Possible choices are {1}'.format(index,', '.join(possible_ninos)))
	lon_name = None
	lat_name = None
	if sst[lon].min() < 0 or sst[lon].max() <= 180:
		lon_name = lon
	if sst[lat][0] > sst[lat][-1]:
		lat_name = lat
	if lon_name is not None or lat_name is not None:
		#print('WARNING: re-arranging SST to be in domain [0,360] x [-90,90]')
		sst = StandardGrid(sst,lon_name,lat_name)

	sst = sst.resample({time:'1M'}).mean()

	def NinoAvg(sst,nino,time,avg,std=True):
		ssta = sst.sel(indList[nino]).mean(dim=[lon,lat])
		sstc = ssta.groupby('.'.join([time,'month'])).mean(dim=time)
		ssta = ssta.groupby('.'.join([time,'month'])) - sstc
		if avg is not None:
			ssta = ssta.load()
			ssta = ssta.rolling({time:avg},min_periods=1).mean()
		if std:
			ssta = ssta/ssta.std(dim=time)
		return ssta

	combined = False
	if isinstance(indList[index],str):
		comb = indList[index]
		if '-' in comb:
			combined = True
		else:
			index = comb
			combined = False

	if combined:
		ind1 = comb.split('-')[0].replace(' ','')
		ind2 = comb.split('-')[1].replace(' ','')
		std = 'std' in ind1
		ind1 = NinoAvg(sst,ind1.replace('std',''),time,None,std)
		std = 'std' in ind2
		ind2 = NinoAvg(sst,ind2.replace('std',''),time,None,std)
		indx = (ind1-ind2).load()
		indx = indx.rolling({time:avg},min_periods=1).mean()
		out  = (indx)/indx.std(dim=time)
		out.name = index.upper()
		return out
	else:
		out = NinoAvg(sst,index,time,avg)
		out.name = index.upper()
		return out


#######################################################
def RollingMeanStd(x,mean_std,r=31,dim='time'):
	'''Compute climatological standard deviation or mean, smoothed by a rolling mean on day of year.

	INPUTS;
	x	: DataArray on which to operate on.
	mean_std: either 'mean' or 'std' depending on what you want to do.
	'''
	import xarray as xr
	import numpy as np
	if mean_std == 'mean':
		xc = x.groupby(dim+'.dayofyear').mean()
	elif mean_std == 'std':
		xc = x.groupby(dim+'.dayofyear').std()
	else:
		raise ValueError("mean_std has to be 'mean' or 'std' but is "+mean_std)
	x1 = xc.roll(dayofyear=2*r,roll_coords=True).rolling(dayofyear=r,center=True).mean().roll(dayofyear=-2*r,roll_coords=True)
	x2 = xc.roll(dayofyear=-2*r,roll_coords=True).rolling(dayofyear=r,center=True).mean().roll(dayofyear=2*r,roll_coords=True)
	xc = xr.where(np.isnan(x1),x2,x1)
	return xc

#######################################################
def Anomaly(da,groupby='time.dayofyear',climate=None):
	'''Compute anomaly of xr.DataArray on frequency given by groupby.

	INPUTS:
	   da	  : xarray.DataArray or xarray.Dataset from which to compute anomalies.
	   groupby: defines frequency for anomalies. If None, return anomalies with respect to total average (useful if only one-dimensional).
	   climate: define period for baseline climatology. list of two dates [string]. All time steps if None.
	OUTPUTS:
	   da	  : input array as anomalies with respect to `groupby`
	'''
	if groupby is None: # only one dimension
		return da - da.mean()
	time_name = groupby.split('.')[0]
	if climate is None:
		clim_filtr = {time_name : slice(da[time_name][0],da[time_name][-1])}
	else:
		clim_filtr = {time_name: slice(str(climate[0]),str(climate[1]))}
	climatology_mean = da.sel(clim_filtr).groupby(groupby).mean(time_name)
	return da.groupby(groupby) - climatology_mean

#######################################################
def Standardize(da,groupby='time.dayofyear',std=None,climate=None):
	'''Standardize xr.DataArray on a frequency given by groupby.
		This is from http://xarray.pydata.org/en/stable/examples/weather-data.html

          INPUTS:
            da     : xarray.DataArray or xarray.Dataset to standardize
            groupby: frequency/basis on which to compute anomalies and standardize. If None, standardize over all dimensions (useful if only one dimensional).
            std    : use this standard deviation to standardize instead of da's standard deviation
 	    climate: use this time period to compute the baseline climatology. Entire time period if None.
	  OUTPUTS:
	    da	   : standarized da
	'''
	if groupby is None: # only one dimension
		return (da - da.mean())/da.std()
	if '.' in groupby:
		time_name = groupby.split('.')[0]
	else:
		time_name = groupby
	if climate is None:
		clim_filtr = {time_name : slice(da[time_name][0],da[time_name][-1])}
	else:
		clim_filtr = {time_name: slice(str(climate[0]),str(climate[1]))}
	climatology_mean = da.sel(clim_filtr).groupby(groupby).mean(time_name)
	if std is None:
		climatology_std = da.sel(clim_filtr).groupby(groupby).std(time_name)
	else:
		climatology_std = std
	stand_anomalies = (da.groupby(groupby) - climatology_mean).groupby(groupby)/climatology_std
	return stand_anomalies
#######################################################
def StackArray(x,dim):
	'''return stacked array with only one dimension left

	INPUTS:
	   x  : xarray.DataArray or Dataset to be stacked
	   dim: sole dimension to remain after stacking, None if all dimensions to be stacked
	OUTPUTS:
	   stacked: same as x, but stacked
	'''
	if len(x.dims) == 0:
		return x
	if dim not in x.dims or dim is None:
		dims = list(x.dims)
	else:
		dims = []
		for d in x.dims:
			if d != dim:
				dims.append(d)
	return x.stack(stacked=dims)

def ComputeStat(i,sx,y,sy,test):
	'''This is part of StatTest, but for parmap.map to work, it has to be an
		independent function.
	'''
	from xarray import DataArray
	from scipy import stats
	# if sx and sy are the same, null hypothesis of the two being distinct
	#  should be rejected at all levels, i.e. p-value of them being
	#  equal should be 100%
	if isinstance(y,DataArray) and sx.shape == sy.shape and (np.all(sx.isel(stacked=i) == sy.isel(stacked=i))):
		ploc = 1.0
	else:
		if test == 'KS':
			_,ploc = stats.ks_2samp(sx.isel(stacked=i),sy.isel(stacked=i))
		elif test == 'MW':
			_,ploc = stats.mannwhitneyu(sx.isel(stacked=i),sy.isel(stacked=i),alternative='two-sided')
		elif test == 'WC':
			_,ploc = stats.wilcoxon(sx.isel(stacked=i),sy.isel(stacked=i))
		elif test == 'T':
			_,ploc = stats.ttest_1samp(sx.isel(stacked=i),y)
		elif test == 'sign': # not really a sig test, just checking sign agreement
                        # note that here a high p-value means significant, as it means
                        #  that a lot of members have the same sign
			lenx = len(sx.isel(stacked=i))
			if y is None: # check sx for same sign
				posx = np.sum(sx.isel(stacked=i) > 0)
				# need to check positive and negative, to account for NaN
				negx = np.sum(sx.isel(stacked=i) < 0)
				ploc = max(posx,lenx-posx)/lenx
			elif isinstance(y,float) or isinstance(y,int) or ( isinstance(sy,DataArray) and 'stacked' not in sy.dims ): # check sx for sign of y
				samex = np.sum( np.sign(sx.isel(stacked=i)) == np.sign(y) )
				ploc  = samex/lenx
			elif isinstance(sy,DataArray) and (len(sy.isel(stacked=i).dims) == 0 ): # check sx for sign of y, where y is not an ensemble, but a function of space
				samex = np.sum( np.sign(sx.isel(stacked=i)) == np.sign(sy.isel(stacked=i)) )
				ploc  = samex/lenx
			else: # check two ensembles for same sign
				# ensembles are not 1-by-1, so we can't check sign along dimension
				# sy might not be an ensemble
				lenx = len(sx.isel(stacked=i))
				posx = np.sum(sx.isel(stacked=i) > 0)/lenx
				leny = len(sy.isel(stacked=i))
				posy = np.sum(sy.isel(stacked=i) > 0)/leny
				ploc = min(posx,posy)/max(posx,posy)
	return ploc

def StatTest(x,y,test,dim=None,parallel=False):
	'''Compute statistical test for significance between
	   two xr.DataArrays. Testing will be done along dimension with name `dim`
	   and the output p-value will have all dimensions except `dim`.

	   INPUTS:
	      x	 : xr.DataArray for testing.
	      y	 : xr.DataArray or scalar for testing against. Or None for single-ensemble sign test.
	      dim: dimension name along which to perform the test.
	      test:which test to use:
		    'KS' -> Kolmogorov-Smirnov
		    'MW' -> Mann-Whitney
		    'WC' -> Wilcoxon
		    'T'  -> T-test 1 sample with y=mean
                    'sign'->test against sign only.
		  parallel: Run the test in parallel? Requires the parmap package.
	   OUTPUTS:
	      pvalx: xr.DataArray containing the p-values.
		     Same dimension as x,y except `dim`.
	'''
	from xarray import DataArray
	if dim is None or len(x.dims) == 1:
		sx = x.expand_dims(stacked=[0])
		parallel = False
	else:
		sx = StackArray(x,dim)
	if parallel:
		import parmap
	nspace = len(sx.stacked)
	if isinstance(y,DataArray):
		if dim is None or len(y.dims) == 1:
			sy = y.expand_dims(stacked=[0])
		else:
			sy = StackArray(y,dim)
	else:
		sy = None
	if parallel:
		pval = parmap.map(ComputeStat,list(range(nspace)),sx,y,sy,test)
	else:
		pval = np.zeros(sx.stacked.shape)
		for i in range(nspace):
			pval[i] = ComputeStat(i,sx,y,sy,test)
	if nspace > 1:
		pvalx = DataArray(pval,coords=[sx.stacked],name='pval').unstack('stacked')
	else:
		pvalx = pval[0]
	return pvalx

#######################################################
def CheckSign(ens,dim,perc):
	'''Check significance using agreement on sign only.
	    This tests ensemble `ens` across dimension `dim` and counts
	    how many members have the same sign. If more than `perc`% of
	    members have the same sign, return is True, else False.

	INPUTS:
	  ens:	the ensemble to test. xr.DataArray or xr.Dataset
	  dim:	dimension to test across. this is normally the ensemble members
	  perc: percentile of significance. require `dim`% of members to have same sign.

	OUTPUTS:
	  filtr: mask of booleans. True means at least `perc`% of members have same sign
	'''
	from xarray import where as xwhere
	from xarray import DataArray
	nmembs = len(ens[dim])
	percentile = 0.01*perc*nmembs
	sign = DataArray(np.sign(ens),coords=ens.coords)
	pos = xwhere(sign > 0,1,0).sum(dim)
	neg = xwhere(sign < 0,1,0).sum(dim)
	pfiltr = pos > percentile # positive sign agrees
	nfiltr = neg > percentile # negative sign agrees
	return pfiltr+nfiltr
#######################################################
def LogPlot(ax,log=True,invert=True,format=True,ylim=None):
	'''Invert y-axis, make it log scale, change ytick formatting.
	'''
	import matplotlib
	if ylim is not None:
		ax.set_ylim(*ylim)
	if log:
		if min(ax.get_ylim()) < 0:
			print('Y-axis extends into negative values. Need ylim input to make it log scale.')
			import sys
			sys.exit()
		ax.set_yscale('log')
	if invert:
		ax.invert_yaxis()
	if format:
		ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

#######################################################
def AddColorbar(fig,axs,cf,shrink=0.95,cbar_args=None):
	'''Add a common colorbar to a FacetPlot.
	    Note that this assumes all panels have been plotted with the same
	    color range and levels.

	  INPUTS:
	     fig   : figure object containing facet plot.
	     axs   : list of all axes in the facet plot.
	     cf	   : image object which contains the values.
		    e.g. cf = plt.contourf(...)
	     shrink: factor to make colorbar smaller.
	     cbar_args: arguments to be passed to pyplot.colorbar
	  OUTPUTS:
	     cbar: colorbar object.
	'''
	if not (isinstance(axs,list) or isinstance(axs,np.ndarray)):
		axs = axs.ravel().to_list()
	if cbar_args is None:
		return fig.colorbar(cf, ax=axs, shrink=shrink)
	else:
		return fig.colorbar(cf, ax=axs, shrink=shrink, **cbar_args)
#######################################################
def AddPanelLabels(axs,loc='lower left',xpos=None,ypos=None,va=None,ha=None,fontsize='large'):
	'''Add a), b), ... labels to each panel within a multipanel figure.

	INPUTS:
	axs : axes object, typically from pyplot.subplots() or a FacetPlot
	loc : semantic location:
	       'lower left': left of and below axes, outside plot
	       'upper left': upper left corner, inside plot
	xpos: manually adjust position in x-direction. units are fraction of axis width. overwrites loc.
	ypos: manually adjust position in y-direction. units are fraction of axis height. overwrites loc.
        va  : manually adjust vertical alignment. overwrites loc.
        ha  : manually adjust horizontal alignment. overwrites loc.
	'''
	from matplotlib.pyplot import Axes
	from numpy import array
	locations = {
		'upper left' : {'xpos':0.01,'ypos':0.99,'va':'top','ha':'left'},
		'lower left' : {'xpos':-0.15,'ypos':-0.15,'va':'bottom','ha':'left'},
		'upper right': {'xpos':0.99,'ypos':0.99,'va':'top','ha':'right'},
		'lower right': {'xpos':1.05,'ypos':-0.15,'va':'bottom','ha':'left'},
		}
	locs = locations[loc]
	if xpos is not None:
	    locs['xpos'] = xpos
	if ypos is not None:
	    locs['ypos'] = ypos
	if va is not None:
	    locs['va'] = va
	if ha is not None:
	    locs['ha'] = ha
	if isinstance(axs,Axes):
	    axs = array(axs)
	from string import ascii_lowercase
	if not isinstance(axs,list):
	    axs = axs.flatten()
	for a,ax in enumerate(axs):
	    ax.text(locs['xpos'],locs['ypos'],ascii_lowercase[a]+')',verticalalignment=locs['va'],horizontalalignment=locs['ha'],transform=ax.transAxes,fontsize=fontsize)

#######################################################
def FindCoordNames(ds):
	'''Find the actual dimension names in xr.Dataset or xr.DataArray
	    which correspond to longitude, latitude, pressure

	   INPUTS:
	   	  ds:        xarray.Dataset or DataArray
	   OUTPUTS:
	      dim_names: Dictionary of dimension names.
		  		      Longitude = ds[dim_names['lon']]
					  Latitude  = ds[dim_names['lat']]
					  Pressure  = ds[dim_names['pres']]
	'''
	odims = list(ds.coords)
	ldims = [d.lower() for d in odims]
	dim_names = {}
	# check for longitude
	for lon in ['longitude','lon','xt_ocean','lon_sub1']:
		if lon in ldims:
			indx = ldims.index(lon)
			dim_names['lon'] = odims[indx]
	for lat in ['latitude','lat','yt_ocean','lat_sub1']:
		if lat in ldims:
			indx = ldims.index(lat)
			dim_names['lat'] = odims[indx]
	for plev in ['level','pres','pfull','lev']:
		if plev in ldims:
			indx = ldims.index(plev)
			dim_names['pres'] = odims[indx]
	return dim_names

#######################################################
def CloseGlobe(ds,lon='infer',copy_all=True):
	'''Close global data cyclicly in longitude to avoid whitespace in stereographic projections.
		This code adapted from code provided by Rishav Goyal.

		INPUTS:
			ds:    xarray.DataArray or xarray.Dataset
			lon:   name of longitude dimension, or 'infer' to automatically find it.
			copy_all: also retain coordinates which are not dimensions of da

		OUTPUS:
			data:  data with 1 additional longitude data point.
	'''
	import xarray as xr
	from cartopy.util import add_cyclic_point as cycpt
	if lon == 'infer':
		dim_names = FindCoordNames(ds)
		lon = dim_names['lon']
	if isinstance(ds,xr.Dataset):
		arrays = list(ds.data_vars.keys())
		is_ds = True
	else:
		arrays = [ds.name]
		is_ds = False
	datas = []
	for array in arrays:
		if is_ds:
			da = ds[array]
		else:
			da = ds
		p1,lon1 = cycpt(da , coord=da[lon])
		coords = []
		for coord in da.dims:
			if coord == lon:
				coords.append((lon,lon1))
			else:
				coords.append((coord,da[coord].data))
		data = xr.DataArray(p1,coords=coords,name=da.name)
		data.attrs = da.attrs
		for coord in da.coords:
			if coord in data.coords:
				data[coord].attrs = da[coord].attrs
			elif copy_all:
				data[coord] = da[coord]
		if is_ds:
			datas.append(data)
	if is_ds:
		return xr.merge(datas)
	else:
		return data

#######################################################
def Regress(x,y,dim):
	'''Compute regression of x onto y along dimension.
	 	Assumes x and y are xarray.DataArrays.

		INPUTS:
			x   : first array, which is regressed onto y
			y   : second array, onto which x is regressed
			dim : dimension along which regression should be performed

		OUTPUTS:
			reg : regression of x onto y, in same units as x
	'''
	from xarray import cov
	return cov(x,y,dim)/y.std(dim)

#######################################################
def DetectEvents(ds,thresh,sep=20,period=1,time='time',kind='auto'):
    '''
    Find events below (kind='min') or above (kind='max') a given threshold for a certain period of time. Events are merged into the earlier if less than a given separation time between them.

    INPUTS:
        ds:     xarray.dataarray used to detect events
        thresh: threshold to define events
        sep:    minimum separation of individual events (from end to start) [number of timesteps]
        period: minimum duration of each event, i.e. #time steps above threshold
        time:   name of time dimension
        kind:   find minimum if 'min', maximum if 'max'.
                if 'auto': minimum if thresh <= 0, maximum elsewhise

    OUTPUTS:
        stats: xarray.Dataset of 
          onset_date:    onset dates of each event
          end_date:      end date of each event
          duration:      duration of each event, i.e. end_date - onset_date
          extreme_value: most extreme value during each event period.
    '''
    import xarray as xr
    from datetime import timedelta
    ds = ds.rename({time:'time'})
    da = ds.rolling(time=period)
    if kind == 'auto':
        if thresh <= 0:
            kind = 'min'
        else:
            kind = 'max'
    if kind.lower() == 'max':
        da = da.min()
        times = da > thresh
    elif kind.lower() == 'min':
        da = da.max()
        times = da < thresh
    # in case no event has been detected
    if np.sum(times) == 0:
        return None
    event_all = da.isel(time=times)
    timestep = ds.time[1] - ds.time[0]
    e = 0
    event_id = [e]
    tm1 = event_all.time[0]
    for t in event_all.time[1:]:
        delta_t = t - tm1
        if delta_t > sep*timestep:
            e += 1
        tm1 = t
        event_id.append(e)
    ex = xr.DataArray(event_id,coords=[event_all.time],name='event_id')
    unique_events = np.unique(ex)
    ts_secs = timedelta(seconds=timestep.values/np.timedelta64(1,'s'))
    zero = 0*ts_secs
    duration = []
    extreme  = []
    start_dates = []
    end_dates = []
    for ue in unique_events:
        filtr = ex == ue
        # this is a workaround to avoid having objects or arrays of arrays
        exdate = ex.isel(time=filtr).time[0].values + zero
        exind  = np.argwhere(ds.time.values==exdate)[0][0]
        prev_date = ds.isel(time=exind-period+1).time.values + zero
        end_date  = ex.isel(time=filtr).time[-1].values + zero
        start_dates.append(prev_date)
        end_dates.append(end_date)
        duration.append(period+sum(filtr.values)-1)
        if kind == 'min':
            extreme.append(ds.sel(time=slice(prev_date,end_date)).min().values)
        elif kind == 'max':
            extreme.append(ds.sel(time=slice(prev_date,end_date)).max().values)
    unicord = [('event',unique_events)]
    onx  = xr.DataArray(start_dates,coords=unicord,name='onset_date')
    edx  = xr.DataArray(end_dates,coords=unicord,name='end_date')
    durx = xr.DataArray(duration,coords=unicord,name='duration')
    extx = xr.DataArray(extreme,coords=unicord,name='extreme_value')
    outx = xr.merge([durx,extx,onx,edx])
    outx.attrs['variable'] = ds.name
    options = {'max':'above','min':'below'}
    outx.attrs['method'] = 'individual {0}-day periods of {4} {1} {2}. Events are considered the same if spaced by less than {3} days'.format(period,options[kind],thresh,sep,ds.name)
    return outx 

#######################################################
def Vorticity(u,v,use_windspharm=False,lon='infer',lat='infer'):
    '''Compute relative vorticity zeta = v_x - u_y on the sphere.
    
     INPUTS:
        u     : zonal wind. xarray.DataArray
        v     : meridional wind. xarray.DataArray
        use_windspharm: whether or not to use windspharm to compute vorticity
        lon   : name of longitude coordinate. Guessed if 'infer'
        lat   : name of latitude coordinate. Guessed if 'infer'

     OUTPUT:
        relative vorticity [1/s]
    '''
    if use_windspharm:
        from windspharm.xarray import VectorWind
        wv = VectorWind(u,v)
        return wv.vorticity()
    else:
        from .constants import a0,coslat
        if lon == 'infer' or lat == 'infer':
            dim_names = FindCoordNames(u)
            if lon == 'infer':
                lon = dim_names['lon']
            if lat == 'infer':
                lat = dim_names['lat']
        one_over_a = 1/a0
        one_over_cosphi = 1/coslat(u[lat])
        dxv = one_over_a*one_over_cosphi*v.differentiate(lon,edge_order=2)
        dyu = one_over_a*                u.differentiate(lat,edge_order=2)
        vort = dxv - dyu
        vort.attrs['units'] = '1/s'
        vort.name = 'relative_vorticity'
        return vort*180/np.pi        

#######################################################
def PotentialTemperature(t,p,p0=1e3):
    '''Compute potential temperature theta = T*(p0/p)**kappa
   
     INPUTS:
        t : temperature in  Kelvins
        p : pressure in same units as p0 [defaults to hPa]
        p0: reference pressure [1000hPa]

     OUTPUTS:
        potential temperature in Kelvins
    '''
    from .constants import kappa
    return t*(p0/p)**kappa

#######################################################
def PotentialVorticity(u,v,t,use_windspharm=False,lon='infer',lat='infer',pres='infer'):
    '''Computes potential vorticity -g*(zeta + f)theta_p. Pressure is assumed in hPa.
    
     INPUTS:
        u     : zonal wind. xarray.DataArray
        v     : meridional wind. xarray.DataArray
        t     : temperature. xarray.DataArray
        use_windspharm: whether or not to use windspharm to compute vorticity
        lon   : name of longitude coordinate. Guessed if 'infer'
        lat   : name of latitude coordinate. Guessed if 'infer'
        pres  : name of pressure coordinate [hPa]. Guessed if 'infer'
  
     OUTPUT:
        potential vorticity. xarray.DataArray
    '''
    from .constants import g,f
    rel_vort = Vorticity(u,v,use_windspharm,lon,lat)
    if lat == 'infer':
        lat = FindCoordNames(u)['lat']
    vort = rel_vort + f(u[lat])
    if pres == 'infer':
        dim_names = FindCoordNames(u)
        pres = dim_names['pres']
    theta = PotentialTemperature(t,t[pres])
    dtheta = theta.differentiate(pres,edge_order=2)*0.01
    return -g*vort*dtheta

#######################################################
def IPV(u,v,t,theta0,use_windspharm=False,lon='infer',lat='infer',pres='infer'):
    ''' Compute Potential Vorticity on Isentropic surface(s).
        Uses wrf-python interplevel function to interpolate from vertical to potential temperature levels.
      
      INPUTS:
        u     : zonal wind. xarray.DataArray
        v     : meridional wind. xarray.DataArray
        t     : temperature. xarray.DataArray
        theta0: value(s) of potential temperature to interpolate onto.
                 this can be a single value or a list/array. 
        use_windspharm: whether or not to use windspharm to compute vorticity
        lon   : name of longitude coordinate. Guessed if 'infer'
        lat   : name of latitude coordinate. Guessed if 'infer'
        pres  : name of vertical coordinate. Guessed (assuming pressure) if 'infer'

      OUTPUT:
	interpolated potential vorticity on given potential temperature levels
    '''
    from wrf import interplevel
    pv = PotentialVorticity(u,v,t,use_windspharm,lon,lat,pres)
    pv = pv.transpose(*u.dims)
    if pres == 'infer':
        pres = FindCoordNames(t)['pres']
    theta = PotentialTemperature(t,t[pres])
    theta = theta.transpose(*u.dims)
    ipv = interplevel(pv,theta,theta0)
    ipv.name = 'ipv'
    ipv.attrs['long_name'] = 'isentropic potential vorticity'
    if 'vert_units' in ipv.attrs:
	    del ipv.attrs['vert_units']
    return ipv


#######################################################
def ComputeSurface(ds,lon='infer',lat='infer'):
    """Compute the surface [m2] of each grid box in ds."""
    from .constants import coslat,a0
    if lon == 'infer':
        lon = FindCoordNames(ds)['lon']
    if lat == 'infer':
        lat = FindCoordNames(ds)['lat']
    #dlam = a0*coslat(ds[lat])*np.deg2rad(np.mean(ds[lon].values[1:]-ds[lon].values[:-1]))
    dlam = a0*coslat(ds[lat])*np.deg2rad(ds[lon].diff(lon))
    #dphi = a0*np.deg2rad(np.mean(ds[lat].values[1:]-ds[lat].values[:-1]))
    dphi = a0*np.deg2rad(ds[lat].diff(lat))
    return dlam*dphi

#######################################################
def GlobalMass(da,lon='infer',lat='infer',pres='infer'):
    """Compute global mass of da, which is assumed in units of [kg/kg], such as specific humidity. Also assumes pressure coordinates in hPa.

	INPUTS:
	  da   - xarray.DataArray assumed in [kg/kg]
	  lon  - name of longitude. Guessed if 'infer'. 
	  lat  - name of latitude. Guessed if 'infer'.
	  pres - name of pressure, assumed in hPa. Guessed if 'infer'.
	OUTPUTS:
	  mass - globally integrated mass
    """
    from .constants import a0,g,coslat
    if lon == 'infer':
        lon = FindCoordNames(da)['lon']
    if lat == 'infer':
        lat = FindCoordNames(da)['lat']
    if pres == 'infer':
        pres = FindCoordNames(da)['pres']
    mass_y = np.deg2rad((coslat(da[lat])*da).integrate(lat))
    mass_x = np.deg2rad(mass_y.integrate(lon))
    mass_p = mass_x.integrate(pres)*100
    mass = a0**2/g*mass_p
    if da[pres][0]>da[pres][-1]:
        mass = -mass
    return mass

#######################################################
def TotalColumnOzone(o3,t,units=1,pres='infer'):
    '''Compute total column ozone in DU, given vertical profile(s) of ozone.
       Ozone concentration is assumed to be mol/mol (units=1) or kg/kg (units=2).
       Assumes pressure in hPa.
    '''
    from .constants import R,Rd,g,Na
    if pres == 'infer':
        pres = FindCoordNames(o3)['pres']
    P = o3[pres]
    if units == 1: #[mol/mol]
        mol_air = P/R/t*100 #[mol/m3]
        mol_o3  = mol_air*o3 #[mol/m3]
        #pressure can be in any direction, so use abs()
        num_o3  = np.abs(Rd/g*(t*mol_o3/P).integrate(pres)) #[mol/m2]
        num_o3  = Na*num_o3 #[molecules/m2]
        one_DU  = 2.69e20 #[molecules/m2]
        DU      = num_o3/one_DU
    return DU
            
def ComputeOzoneHoleArea(o3,hemi='S',lat='infer',lon='infer'):
        if lat == 'infer':
                lat = FindCoordNames(o3)['lat']
        if lon == 'infer':
                lon = FindCoordNames(o3)['lon']
        o3 = StandardGrid(o3)
        if hemi == 'S':
                o3s = o3.sel({lat:slice(-90,-40)})
        else:
                o3s = o3.sel({lat:slice(40,90)})
        filtr = o3s < 220
        surf = ComputeSurface(o3s)
        return surf.where(filtr).sum(['lon','lat'])
