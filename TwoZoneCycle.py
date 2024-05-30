#%% LIBRARIES
import numpy as np
import matplotlib.pyplot as plt
import scipy
from Thermodynamic import Combustion, Equilibrium
from BurnedFraction import x_Wiebe
from Volume import volume
from Twall import Tw
from Woschni import hc
from AdiabaticTemperature import T_ad
from Zeldovich import zeldovich

def homogeneous (Variables):  
    def integrate( THETA, THETAE, Y):
        def rates( THETA, Y ):
            YPRIME = np.zeros(NY)
            VOL = volume(THETA)[0]
            X = x_Wiebe(THETA,THETAS,THETAB, PHI) [0]
            EM = np.exp(-BLOWBY*(THETA*np.pi/180 + np.pi)/OMEGA)
            M = EM*MNOT
            VOL, DV, Ac = volume(THETA)
            AA = (DV + VOL*BLOWBY/OMEGA)/M
            C1 = HEAT*(Ac)/OMEGA/M/1000
            C0 = np.sqrt(X)
            TW = Tw(THETA)
            
            
            P = Y[0]
            TB = Y[1]
            TU = Y[2]
            
            if ( X > 0.999 ):
                ierr, YB = Equilibrium( TB, P, PHI)
                if ( ierr != 0 ):
                    print('Error in Equilibrium({0}, {1}, {2}): {3}\n'.format(TB, P, PHI, ierr))
                xxx, HL, xxx, xxx, VB, xxx, CP, xxx, DVDT, DVDP = Combustion( TB, P, PHI, f)
                     
                
                BB = C1/CP*DVDT*TB*(1-TW/TB)
                CC = 0
                DD = 1/CP*TB*DVDT**2 + DVDP
                EE = 0
                
                YPRIME[0] = (AA + BB + CC)/(DD + EE)
                YPRIME[1] = -C1/CP*(TB-TW) + 1/CP*DVDT*TB*YPRIME[0]
                YPRIME[2] = 0
                
            elif ( X > 0.001 ):
                ierr, YB = Equilibrium( TB, P, PHI)
                if ( ierr != 0 ):
                    print('Error in Equilibrium({0}, {1}, {2}): {3}\n'.format(TB, P, PHI, ierr))
                xxx, HU, xxx, xxx, VU, xxx, CPU, xxx, DVDTU, DVDPU = Combustion( TU, P, PHI, f)
                xxx, HB, xxx, xxx, VB, xxx, CPB, xxx, DVDTB, DVDPB = Combustion( TB, P, PHI, f)
                  
                BB = C1*(1/CPB*TB*DVDTB*C0*(1-TW/TB) + 1/CPU*TU*DVDTU*(1-C0)*(1-TW/TU))
                DX = x_Wiebe(THETA,THETAS,THETAB, PHI)[1]
                CC = -(VB-VU)*DX - DVDTB*(HU-HB)/CPB*(DX-(X-X**2)*BLOWBY/OMEGA)
                DD = X*(VB**2/CPB/TB*(TB/VB*DVDTB)**2 + DVDPB)
                EE = (1-X)*(1/CPU*TU*DVDTU**2 + DVDPU)
                HL = (1-X**2)*HU + X**2*HB
                
                YPRIME[0] = (AA + BB + CC)/(DD + EE)
                YPRIME[1] = -C1/CPB/C0*(TB-TW) + 1/CPB*TB*DVDTB*YPRIME[0] + (HU-HB)/CPB*(DX/X - (1-X)*BLOWBY/OMEGA)
                YPRIME[2] = -C1/CPU/(1+C0)*(TU-TW) + 1/CPU*TU*DVDTU*YPRIME[0]
                
            else:           
                xxx, HL, xxx, xxx, xxx, xxx, CP, xxx, DVDT, DVDP = Combustion( TU, P, PHI, f)
                
                BB = C1*1/CP*TU*DVDT*(1-TW/TU)
                CC = 0
                DD = 0
                EE = 1/CP*TU*DVDT**2 + DVDP
                
                YPRIME[0] = ( AA + BB + CC )/(DD + EE)
                YPRIME[1] = 0
                YPRIME[2] = -C1/CP*(TU-TW) + 1/CP*TU*DVDT*YPRIME[0]
            
            YPRIME[3] = Y[0]*DV;
            
            if ( X > 0.001 ):
                for k in range(NNOX) : 
                    if ( THETA >= THETAS + (k-1)/(NNOX-1)*THETAB ):
                        YPRIME[Npara+k] = zeldovich( TB, P/100, YB, Y[Npara+k]/(MW_NO*VB*1000) )*MW_NO*VB*1000/OMEGA
            
            for JJ in range(NY) : 
                YPRIME[JJ] = YPRIME[JJ]*np.pi/180
            
            return YPRIME
        sol = scipy.integrate.solve_ivp(rates, np.array([THETA, THETAE]), Y, method='RK23', first_step = (THETAE-THETA)/100)
        
        return sol
    #%% Variables to make vary
    PHI, P1, T1, RPM, BLOWBY, THETAS, THETAB, f = Variables 

    #%% Global variables
    HEAT = 500
    OMEGA = RPM*np.pi/30

    THETA = -270
    DTHETA = 1
    THETAE = THETA+DTHETA

    VOL = volume(THETA)[0]
    EM = np.exp(-BLOWBY*(THETA*np.pi/180 + np.pi)/OMEGA)
    
    to_ppm = 10**6
    MW_NO = 30
    NNOX = int(THETAB/DTHETA)
    
    Npara = 4
    NY = Npara+int(NNOX)
    Y = np.zeros(NY)
    Y[0] = P1
    Y[1] = 0
    Y[2] = T1

    vU = Combustion(Y[2], Y[0], PHI, f)[4]

    MNOT = VOL/vU

    NN = 540;
    SAVE_THETA = np.zeros(NN)
    SAVE_VOL = np.zeros(NN)
    SAVE_TB = np.zeros(NN)
    SAVE_TU = np.zeros(NN)
    SAVE_P = np.zeros(NN)
    SAVE_X = np.zeros(NN)

    for II in range(NN) :
        SAVE_VOL[II] = volume(THETA)[0]
        SAVE_X[II] = x_Wiebe(THETA,THETAS,THETAB, PHI) [0]  
        EM = np.exp(-BLOWBY*(THETA*np.pi/180 + np.pi)/OMEGA)
        sol = integrate( THETA, THETAE, Y )
        
        for i in range(len(Y)):
            Y[i] = sol.y[i][-1]
        
        if sol.status != 0 :
            print('Integration error fail to reach THETAE = {0:.2f}'.format(THETAE))
            print(SAVE_P,SAVE_TU,SAVE_TB)
            break  
        THETA=THETAE;
        THETAE=THETA+DTHETA;
        
        # save data for plotting late
        SAVE_THETA[II] = sol.t[-1];
        SAVE_P[II] = sol.y[0][-1];
        SAVE_TB[II] = sol.y[1][-1];
        SAVE_TU[II] = sol.y[2][-1];

        if ( THETAS >= THETA) and (THETAS < THETAE ):
            Y[1] = T_ad( Y[0], Y[2], PHI, f)
        
        if ( THETA > THETAS + THETAB ):
            Y[2] = 0
    
    NOX_ppm = 0;
    for nn in range(NNOX) :
        THETA = THETAS + (nn-1)/(NNOX-1)*THETAB
        dxbdtheta = 0.5*np.sin(np.pi*(THETA-THETAS)/THETAB)*np.pi/THETAB
        dxb = dxbdtheta*DTHETA
        NOX_ppm = NOX_ppm + Y[Npara+nn]*dxb*to_ppm
                
    IMEP = Y[3]/(SAVE_VOL[0])
    
    return SAVE_THETA, SAVE_VOL, SAVE_P, SAVE_TB, SAVE_TU, SAVE_X, IMEP, NOX_ppm