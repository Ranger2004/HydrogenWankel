#%% LIBRARY
import numpy as np
R = 8.314

#%% FONCTION
def fuel (T) :
    """
    Computes the properties of hydrogen fuel

    Parameters
    ----------
    T : float
        Temperature - [K]

    Returns
    -------
    alpha : float
        Number of carbon molecules
    beta : float
        Number of hydrogen molecules
    gamma : float
        Number of oxygen molecules
    delta : float
        Number of sulfur molecules
    h : float
        Enthalpy - [kJ/kg]
    s : float
        Entropy - [kJ/kg/K]
    cp : float
        Heat capacity of constant pressure - [kJ/kg/K]
    mw : float
        Molecular weight - [kg/kmol]
    Fs : float
        Air-fuel ratio - [/]
    q : float
        Lower heating value - [kJ/kg]

    """
    if T>1000:
        FuelProps = [.3189e+1  ,.133822e-2,-.528993e-6,.959193e-10,-.648479e-14,
                     .982832e+4,.674581e+1 ]
    
    else :
        FuelProps = [ 0.30574451e+1,   0.26765200e-2,    -0.5099162e-5,  0.55210391e-8,
                 -0.18122739e-11, -0.98890474e+3, -0.22997056e+1]
    # Fuel chemical formula
    #            C         H        O        N
    #            alpha     beta     gamma    delta
    FuelInfo = [ 0,         2,        0,        0 ]
    
    # stoichiometric fuel-air ratio
    FSv = 1/34.33
    
    # available energy of combustion ac
    ac = 120000
    
    # stoichiometric fuel-air ratio
    Fs = FSv
    
    # available energy
    q = ac
    
    # Fuel composition
    alpha = FuelInfo[0]
    beta = FuelInfo[1]
    gamma = FuelInfo[2]
    delta = FuelInfo[3]
    
    # Compute fuel properties
    ao = FuelProps[0]
    bo = FuelProps[1]
    co = FuelProps[2]
    do = FuelProps[3]
    eo = FuelProps[4] 
    
    # Compute thermodynamic properties
    h = ao + bo/2*T +co/3*T**2 +do/T
    s = ao*np.log(T) + bo*T +co/2*T**2 + eo
    cp = ao + bo*T + co*T**2
    
    # Calculate molecular weight of fuel
    mw = 12.01*alpha + 1.008*beta + 16.00*gamma + 14.01*delta
    
    q = q*mw*10**(-3)
    
    return alpha, beta, gamma, delta, h, s, cp, mw, Fs, q