import numpy as np

def Tw(x):
    """
    Conversion of in-wall data into a regressed value for any angle

    Parameters
    ----------
    x : float
        Crank angle - [deg]

    Returns
    -------
    sol : float
        Temperature wall - [K]

    """
    # For 4000 rpm
    T_w = [120,110,100,101,110,130,140,190,160,130,120,120,120]
    for i in range(len(T_w)):
        T_w[i]+=273.15
    theta = [-180,-145,-103,-90,-65,-24,0,35,58,90,120,139,180]
    for i in range(len(theta)):
        theta[i]*=3
    degree = 7
    
    a = np.polyfit(theta,T_w,degree)
    sol = 0
    for i in range(len(a)):
        sol += a[i]*x**(degree-i)
    
    return sol