#%% LIBRARY
import numpy as np
import scipy
from Qin import Q

def HRE (Variables):
    def volume(theta):
       b = 0.080
       e = 0.015
       R = 0.105
       Rx = (R**2-2*e*R+4*e**2)/(R-4*e)
       beta = (np.sqrt(3*R))/((6*e*R)/(R-4*e)+R+2*e)
       L = 2*Rx*beta
       K = R/e
       phi_max = np.arcsin(3/K)
       
       F = np.pi/3*e**2 + e*R*(2*np.cos(phi_max*np.pi/180)-(3*np.sqrt(3))/2*np.sin(2/3*(theta)+np.pi/2)) + (2/9*R**2+4*e**2)*phi_max
       V = b*F
       
       DV = -np.sqrt(3)*e*R*b*np.cos(2/3*theta+np.pi/2)
       
       A = V/b + V/L + b*L
        
       return V, DV, A
    
    def rates(theta,fy) :
        vol, dvol, A = volume(theta)
        dx=0.;
        if(theta>thetas) and (theta<thetas+thetad) :
            dx = 0.5*np.sin(np.pi*(theta-thetas)/thetad )*180/(thetad*180/np.pi)
        
        term1= -gamma*fy[0]*dvol/vol
        term2= (gamma-1)/vol*Qin*dx*10**(-6)
        yprime= term1 + term2
        return yprime
    
    def integrate_ht(theta,thetae,fy) :
        y = scipy.integrate.solve_ivp(rates, np.array([theta, thetae]), fy, method='RK23')
        fy = y;
        return fy
    
    thetas = Variables[3]*np.pi/180
    thetad = Variables[4]*np.pi/180
    gamma = 1.4
    phi = Variables[0]
    Qin = Q(2400,phi)[0]
    T1 = Variables[2]
    
    step=1*np.pi/180
    NN=int(3*np.pi/step)
    
    theta = -270*np.pi/180
    thetae = theta + step
    
    th=np.zeros(NN)
    volu=np.zeros(NN) 
    press=np.zeros(NN)
    temp = np.zeros(NN)
    fy=np.zeros(1)
    fy[0] = Variables[1]
    
    for i in range(NN):
        sol = integrate_ht(theta,thetae,fy)
        volu[i] = volume(theta)[0]
        fy[0] = sol.y[0][-1]
        theta = thetae;
        thetae = theta+step;
        
        th[i]=theta
        press[i]=sol.y[0][-1]
        
    def Temp(p):
        temp = np.zeros(NN)
        temp[0]=T1
        for i in range(1,len(p)):
            temp[i] = temp[i-1]*(p[i-1]/p[i])**((1-gamma)/gamma)
        return temp
    temp = Temp(press)
    return th, press, temp, volu
