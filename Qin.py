#%% LIBRARY
import numpy as np
import pyromat as pm
from Fuel import fuel
H2O = pm.get('ig.H2O')
N2 = pm.get('ig.N2')
O2 = pm.get('ig.O2') 
H2 = pm.get('ig.H2')
CO2 = pm.get('ig.CO2')
CO = pm.get('ig.CO')
sp = [CO2,H2O,N2,O2,CO,H2]
R = 8.314


#%% HEAT OF COMBUSTION CALCULATIONS
def Q(T,phi):
    """
    Computation of the heat released by the combustion of hydrogen

    Parameters
    ----------
    T : float
        Temperature - [K]
    phi : float
        Equivalence ratio - [/]

    Returns
    -------
    qc : float
        Released energy - [kJ]
    m_fa : float
        Fuel-air mixture mass - [kg]

    """
    # On veut la composition du fuel
    alpha, beta, gamma, delta, h_fuel, so_fuel, cp_fuel, m_fuel, Fs, q = fuel(T)
    
    # On veut les masse moleculaires des constituants
    #      CO2     H2O      N2        O2       CO       H2
    Mi = [CO2.mw(),H2O.mw(),N2.mw(),O2.mw(),CO.mw(),H2.mw()]
    
    # Conversion de phi en lambda
    l = 1/phi
    
    # On veut les moles des reactifs
    n_fuel = 1
    n_O2 = l/2
    n_N2 = l/2*3.76
    
    # On veut la masse de fuel -> [mol]*[kg/kmol] = [g]
    m_fa = n_fuel*m_fuel + n_O2*O2.mw() + n_N2*N2.mw()
    
    # On considère qu'on a pas de résidus
    n = np.array([0,(l+1)/2,3.76*l/2,(l-1)/4,0,(1-l)/2])
    m = np.zeros(len(n))
    for i in range(len(n)):
        if n[i]<=0:
            n[i]=0
        m[i] = n[i]*Mi[i]*10**(-3)
    
    qc = n_fuel*Mi[5]*h_fuel + n_O2*Mi[3]*O2.h(T) + n_N2*Mi[2]*N2.h(T)
    
    for i in range(len(n)):
        qc-=m[i]*sp[i].h(T)
    return qc[0]*m_fa*10**(-3), m_fa