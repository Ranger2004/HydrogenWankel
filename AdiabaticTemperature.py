#%% LIBRARY
from Thermodynamic import Combustion

Kelvin = 273.15

#utiliser kwarg ou arg dans argument pour randre condition optionnel

def T_ad( P, TU, PHI, f):
    """
    Iteration process to compute the adiabatic temperature of a mixture proposed by Ferguson.

    Parameters
    ----------
    P : float
        Pressure - [kPa]
    TU : float
        Initial temperature - [K]
    PHI : float
        Equivalence ratio - [/]
    f : float
        Residual factor - [/]
    cond : integer
        Condition of study ;
            0 = general combustion
            1 = simplified combustion        
    Returns
    -------
    TB : float
        Adiabatic flame temperature at constant pressure - [K]

    """
    TB = 2000;
    HU = Combustion( TU, P, PHI, f )[1];
    for ITER in range(50) :
        x, HB,x,x,x,x, CP,x,x,x = Combustion( TB, P, PHI, f );
        DELT = +(HU-HB)/CP;
        TB = TB + DELT;
        
        if ( abs(DELT/TB) < 0.001 ):
            break;
    return TB

