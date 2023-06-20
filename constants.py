g = 9.81 #[m/s2 ]
a0 = 6376.0e3 # [m]
Omega = 7.292e-5 #[1/s]
R  = 8.314 #[J/mol/K]
Rd = 287.04 # [J/kg/K]
Rv = 461.50 # [J/kg/K]
Na = 6.022e23 # [1/mol]
kappa = 2./7. # []
cp = Rd/kappa # [J/kg/K]
sigma = 5.6734e-8 # [W/m2/K4]
ES0 = 610.78 # [Pa]
HLV = 2.5e6 # [J/kg]
Tfreeze = 273.16 # [K]
p0 = 1e3 # [hPa]
p0_Pa = 1e5 # [Pa]

cmaps = {
    'slp'    : 'BrBG_r',
    'precip' : 'PuOr',
    't'      : 'RdBu_r',
    'u'      : 'RdBu_r',
    'olr'    : 'PiYG',
    'cldfrac': 'binary_r',
}

def f(lat):
    ''' Compute Coriolis parameter f = 2*Omega*sin(lat)

    INPUTS:
      lat:  latitude in degrees. Either numpy or xarray array
    OUTPUTS:
      f:    Coriolis parameter. Either numpy or xarray array
    '''
    from numpy import sin,deg2rad
    return 2*Omega*sin(deg2rad(lat))

def beta(lat,u=None):
    ''' Compute meridional derivative of Coriolis parameter, i.e. beta = 2*Omega*cos(lat)/a0.
          If u is not None, returns beta = 2*Omega*cos(lat)/a0 - u_yy
          If u is not None, it has to be an xarray.DataArray.

    INPUTS:
      lat:  latitude in degrees. Either numpy or xarray array
      u:    None or meridonal wind as xarray dataarray. If
    OUTPUTS:
      beta:   beta parameter. Either numpy or xarray array depending on input.
    '''
    from numpy import cos,deg2rad
    beta = 2*Omega*cos(deg2rad(lat))/a0
    if u is None:
        return beta
    else:
        lats = lat.name
        # uyy = d_phi(1/acosphi*d_phi(u*cosphi))/a0
        uy = deg2rad((u*coslat(lat))).differentiate(lats,edge_order=2)
        uyy = deg2rad((uy/coslat(lat)/a0).differentiate(lats,edge_order=2))/a0
        return beta - uyy


def coslat(lat):
    '''Compute cosine of latitude from degrees.
    '''
    from numpy import cos,deg2rad
    return cos(deg2rad(lat))

def sinlat(lat):
    '''Compute sine of latitude from degrees.
    '''
    from numpy import sin,deg2rad
    return sin(deg2rad(lat))
