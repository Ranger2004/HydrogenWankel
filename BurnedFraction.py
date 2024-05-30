import numpy as np

#%% FUNCTION

def xb(THETA,THETAS,THETAB):
    """
    Computes the burned fraction and its derivative.

    Parameters
    ----------
    THETA : float
        Crank angle - [deg]
    THETAS : float
        Start of combustion angle - [deg]
    THETAB : float
        Duration of combustion angle - [deg]

    Returns
    -------
    X : float
        Burned fraction
    DX : float
        Derivative of the burned fraction according to the crank angle

    """
    X = 0.5*(1-np.cos(np.pi*(THETA-THETAS)/THETAB));
    DX = 0.5*np.sin(np.pi*(THETA-THETAS)/THETAB)*180/(THETAB)
    if ( THETA <= THETAS ):
        X = 0.
        DX = 0
    if ( THETA >= THETAS+THETAB ):
        X = 1.
        DX = 0
    return X,DX

def Wiebe_m(PHI):
    """
    Regression for the Wiebe form factor

    Parameters
    ----------
    PHI : float
        Equivalence ratio

    Returns
    -------
    sol : float
        Wiebe's form factor

    """
    n = [1.05,1.07,1.1,1.2]
    phi = [1/1.4, 1/1.6, 1/1.8, 1/2]
    degree = 1
    
    a = np.polyfit(phi,n,degree)
    sol = 0
    for i in range(len(a)):
        sol += a[i]*PHI**(degree-i)
    return sol


def x_Wiebe(THETA,THETAS,THETAB, PHI):
    """
    Compute the burned fraction rate according to Wiebe's formula

    Parameters
    ----------
    THETA : float
        Crank angle - [deg]
    THETAS : float
        Start of combustion angle - [deg]
    THETAB : float
        Duration of combustion angle - [deg]
    PHI : float
        Equivalence ratio

    Returns
    -------
    X : float
        Burned fraction
    DX : float
        Derivative of the burned fraction according to the crank angle

    """
    a = 4
    m = Wiebe_m(PHI)
    X = 1-np.exp(-a*((THETA-THETAS)/THETAB)**(m+1));
    DX = (m+1)*a*(1-X)*((THETA-THETAS)/THETAB)**(m)
    if ( THETA <= THETAS ):
        X = 0.
    if ( THETA >= THETAS+THETAB ):
        X = 1.
    return X,DX
    