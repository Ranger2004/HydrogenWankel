import numpy as np

def zeldovich(T, p, y, NO):
    """
    Computing the Zeldovich model extended by Lavoie, Heywood and Keck in 1970

    Parameters
    ----------
    T : float
        Temperature - [K]
    p : float
        Pressure - [bar]
    y : float
        Specific volume of the mixture
    NO : float
        Concentration of nitric oxide already present

    Returns
    -------
    dNOdt : float
        The rate of change of nitric oxide

    """
    k1 = 1.8*10**14*np.exp(-38370/T);
    k2r = 3.8*10**9*T*np.exp(-20820/T);
    k3r = 1.7*10**14*np.exp(-24560/T);
    
    N_V = (100000*p)/(8.314*T)*(1/100)**3;
    
    N2e = y[1]*N_V;
    He = y[4]*N_V;
    Oe = y[5]*N_V;
    NOe = y[6]*N_V;
    
    R1 = k1*Oe*N2e;
    R2 = k2r*NOe*Oe;
    R3 = k3r*NOe*He;

    alpha = NO/NOe;
    dNOdt = 2*R1*(1-alpha*alpha)/(1+alpha*R1/(R2+R3));
    return dNOdt
