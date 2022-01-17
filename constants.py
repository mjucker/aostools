g = 9.81 #[m/s2 ]
a0 = 6376.0e3 # [m]
Omega = 7.292e-5 #[1/s]
Rd = 287.04 # [J/kg/K]
Rv = 461.50 # [J/kg/K]
kappa = 2./7. # []
cp = Rd/kappa # [J/kg/K]
sigma = 5.6734e-8 # [W/m2/K4]
ES0 = 610.78 # [Pa]
HLV = 2.5e6 # [J/kg]
Tfreeze = 273.16 # [K]
p0 = 1e3 # [hPa]
p0_Pa = 1e5 # [Pa]

def f(lat):
    from numpy import sin,deg2rad
    return 2*Omega*sin(deg2rad(lat))

def coslat(lat):
    from numpy import cos,deg2rad
    return cos(deg2rad(lat))
