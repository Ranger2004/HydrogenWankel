import numpy as np
from FiniteEnergyRelease import HRE
from BurnedFraction import xb

#Engine data
b = 80
e = 15
Rad = 105
Rx = (Rad**2-2*e*Rad+4*e**2)/(Rad-4*e)
beta = (np.sqrt(3*Rad))/((6*e*Rad)/(Rad-4*e)+Rad+2*e)
L = 2*Rx*beta
K = Rad/e
phi_max = np.arcsin(3/K)
OMEGA = 4000*np.pi/30


def hc (Variables):    
    thetas = Variables[3]
    thetad = Variables[4]
    th,P,T,vol = HRE(Variables)[:4]
    index = 0
    x = np.zeros(len(th))
    for i in range(len(th)):
        x[i] = xb(th[i]*180/np.pi, thetas, thetad)[0]
        if x[i]>0:
            index = i-1
    V1 = vol[index]
    p1 = P[index]
    T1 = T[index]
    pm = 20
    
    m = 0.8
    C = 0.013
    C1 = 2.28
    C2 = 0.00324
    
    Vs = max(vol)
    
    V = vol
    d = 2/(1/b + b*L/V)
    cm = Rx*OMEGA/3
    wc = C2*Vs*T1/(p1*V1)*(P-pm)
    w = C1*cm+wc
    
    return C*d**(m-1)*(P)**m*w**m*T**(0.75-1.62*m)
