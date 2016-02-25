#!/usr/bin/python
# Filename: test_climate
import pytest

# get the functions from pv_tools.climate
execfile('../climate.py')

@pytest.fixture
def testFile():
    import os
    workDir = os.getcwd()+'/'
    return workDir,'climate_test.nc'
#######################################################


# now define the tests
def test_CheckAny():
    # one word within long string
    assert CheckAny('abc','abc def ghi') == 1
    # one word of a set of words within long string
    assert CheckAny(('aaa','def','ggg'),'abc def ghi') == 1
    # word not contained in long string
    assert CheckAny(('aaa','deg','aeh','efg','iab'),'abc def ghi') == 0

def test_FindDayOfYear(testFile):
    import os
    workDir = testFile[0]
    inFile = testFile[1]
    calendar_types = ['standard', 'gregorian', 'proleptic_gregorian', 'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day']
    monthList=['Ja','Fe','Ma','Ap','My','Jn','Jl','Au','Se','Oc','No','De','DJF','JJA','annual']
    for calendar in calendar_types:
        for month in monthList:
            outFile = ComputeClimate(inFile,climatType=month,wkdir=workDir,cal=calendar)
            assert os.path.isfile(outFile)
            os.remove(outFile)

def test_ComputeSaturationMixingRatio(testFile):
    import netCDF4 as nc
    inFile = nc.Dataset(''.join(testFile),'r')
    T = inFile.variables['temp'][:]
    p = inFile.variables['pfull'][:]
    qsat = ComputeSaturationMixingRatio(T,p,1)
    assert T.shape == qsat.shape
    # check a random point
    import numpy.random as rd
    from numpy import exp
    i = rd.random_integers(T.shape[0]-1)
    j = rd.random_integers(T.shape[1]-1)
    k = rd.random_integers(T.shape[2]-1)
    t = rd.random_integers(T.shape[3]-1)
    esat = 610.78*exp(2.5e6*(1./273.16 - 1/T[i,j,k,t])/461.5)
    esat = 287.04/461.5*esat/(p[j]*100-esat)
    assert abs(qsat[i,j,k,t]-esat) < 1.e-9

def test_ComputePsi(testFile):
    def checkPsi(data,outFile,shpe):
        psi,psis = ComputePsi(data,outFile)
        assert psi.shape  == shpe[:-1]
        assert psis.shape == shpe[:-1]
        if outFile is not 'none':
            import os
            assert os.path.isfile(outFile)
            os.remove(outFile)
    # first, test ComputePsi with input variables
    import netCDF4 as nc
    fileName = ''.join(testFile)
    inFile = nc.Dataset(fileName,'r')
    data = {}
    data['temp'] = inFile.variables['temp'][:]
    data['vcomp']= inFile.variables['vcomp'][:]
    data['lat']  = inFile.variables['lat'][:]
    data['pfull']= inFile.variables['pfull'][:]
    shpe = data['temp'].shape
    checkPsi(data,'none',shpe)
    #next, test with input file
    checkPsi(fileName,testFile[0]+'test.nc',shpe)

def test_ComputeVertEddy(testFile):
    def CheckVertEddy(v,t,p,w):
        shpe = v.shape
        v_bar,t_bar = ComputeVertEddy(v,t,p,p0=1e3,wave=w)
        assert v_bar.shape == shpe[:-1]
        assert t_bar.shape == shpe[:-1]
    import netCDF4 as nc
    inFile = nc.Dataset(''.join(testFile),'r')
    v = inFile.variables['vcomp'][:]
    t = inFile.variables['temp'][:]
    p = inFile.variables['pfull'][:]
    nModes = v.shape[-1]/2
    # check total value
    CheckVertEddy(v,t,p,-1)
    from numpy.random import random_integers
    w = random_integers(nModes)-1
    CheckVertEddy(v,t,p,w)

def test_eof(testFile):
    assert True
