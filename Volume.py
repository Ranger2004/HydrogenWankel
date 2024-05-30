import numpy as np

def volume(theta) :  
    """
    Computation of the combustion chamber characteristics

    Parameters
    ----------
    theta : float
        Crank angle - [deg]

    Returns
    -------
    V : float
        Volume - [$m^3$]
    DV : float
        Dv/dtheta - [$m^3$/rad]
    Ac : float
        In-wall area - [$m^2$]

    """
    theta = theta*np.pi/180
    b = 0.080
    e = 0.015
    R = 0.105 - 0.01 # Correction factor
    Rx = (R**2-2*e*R+4*e**2)/(R-4*e)
    beta = (np.sqrt(3*R))/((6*e*R)/(R-4*e)+R+2*e)
    L = 2*Rx*beta
    K = R/e
    phi_max = np.arcsin(3/K)
    
    F = np.pi/3*e**2 + e*R*(2*np.cos(phi_max*np.pi/180)-(3*np.sqrt(3))/2*np.sin(2/3*(theta)+np.pi/2)) + (2/9*R**2+4*e**2)*phi_max
    V = b*F
    
    DV = (-np.sqrt(3)*e*R*b*np.cos(2/3*theta+np.pi/2))
    
    Ac = V/b + V/L + b*L
    
    return V, DV, Ac