#%% LIBRARY
import numpy as np
from Fuel import fuel

#%% FONCTION
def Combustion( T, p, phi, f):
    """
    Computes the thermodynamic properties for the simplified combustion equation,
    H2 + l/2*(O2 + 3.76*N2) -> (l+1)/2*H2O + (1-l)/2*H2 + (l-1)/4*O2 + 3.76*l/2*N2

    Parameters
    ----------
    T : integer
        Temperature - [K].
    p : integer
        Pressure - [kPa].
    phi : integer
        Richness of the fresh gas mixture - [/].
    f : integer
        Residual fraction - [/].

    Returns
    -------
    Y : array
        Mole fraction of constituents
            y(0)  : H2O
            y(1)  : N2
            y(2)  : O2
            y(3)  : H2
    h : integer
        Specific enthalpy - [kJ/kg].
    u : integer
        Specific intenal energy - [kJ/kg].
    s : integer
        Specific entropy - [kJ/kg/K].
    v : integer
        Specific volume - [m3/kg].
    R : integer
        Specific ideal gas constant - [kJ/kg/K].
    Cp : integer
        Specific heat capacity at constant pressure - [kJ/kg/K].
    MW : integer
        Molecular weight of the mixture - [kg/kmol].
    dvdT : integer
        (dv/dT) at const P - [m3/kg/K].
    dvdP : TYPE
        (dv/dP) at const T - [m3/kg/kPa].

    """
    
    # Composition of fuel
    alpha, beta, gamma, delta, h_fuel, so_fuel, cp_fuel, m_fuel, Fs, q = fuel(T)
    
    # Curve fit coefficients for thermodynamic properties
    # Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    # h/RT = a1 + a2/2*T + a3/3*T^2 + a4/4*T^3 + a5/5*T^4 + a6/T
    # so/R = a1*ln(T) + a2*T + a3/2*T^2 + a4/3*T^3 + a5/4*T^4 + a7
    if ( T >= 1000 ) :
        # high temp curve fit coefficients for thermodynamic properties 1000 < T < 3000 K
        A = [[2.67703787e+0,2.97318329e-3,-7.7376969e-7,9.44336689e-11, -4.26900959e-15,-2.98858938e+4,6.88255571e+0 ],      # H2O
             [2.95257626e+0,1.39690057e-3,-4.92631691e-7,7.86010367e-11,-4.60755321e-15,-9.23948645e+2,5.87189252e+0 ],      # N2
             [3.66096083e+0,6.56365523e-4,-1.41149485e-7,2.05797658e-11,-1.2991324e-15,-1.21597725e+3,3.41536184e+0 ],       # O2
             [2.93286579e+0,8.26607967e-4,-1.46402335e-7,1.54100359e-11,-6.88804432e-16,-8.13065597e+2,-1.02432887e+0 ] ]    # H2
    
    if T < 1000 :
        # low temp curve fit coefficients for thermodynamic properties, 300 < T <= 1000 K
        A = [[ 4.19864056e+0, -2.03643410e-3,  6.52040211e-6, -5.48797062e-9, 1.77197817e-12, -3.02937267e+4, -8.49032208e-1 ], # H2O
             [ 3.53100528e+0, -1.23660987e-4,  -5.02999437e-7, 2.435306e-9, -1.40881235e-12, -1.04697628e+3,  2.96747468e+0 ], # N2
             [ 3.78245636e+0, -2.99673415e-3,  9.84730200e-6, -9.68129508e-9, 3.24372836e-12, -1.06394356e+3, 3.65767573e+0 ], # O2
             [ 2.34433112e+0, 7.98052075e-3, -1.94781510e-5,  2.01572094e-8, -7.37611761e-12, -9.17935173e+2, 6.83010238e-1 ] ] # H2
        
    # molecular weights of constituents (g/mol)
    #      H2O      N2         O2      H2
    Mi = [18.02,   28.013,   32.00,  2.016 ];
    
    # air-fuel ratio stoichiometric
    a_s = alpha + beta/4 - gamma/2;
    
    # mole fraction of fuel, O2, N2
    y_1 = 1 / (1 + 4.76*a_s/phi);  # mole fraction for one mole of reactant
    y_fuel = y_1; # assuming 1 mole fuel
    y_O2 = a_s/phi * y_1; # a_s/phi moles O2
    y_N2 = a_s/phi*3.76 * y_1; # a_s/phi * 3.76 moles N2
    
    # On veut la masse de fuel
    m_fa = y_fuel*m_fuel + y_O2*32.00 + y_N2*28.013;
    
    Y = np.zeros(4) #[H2O,N2,O2,H2]
    n = np.zeros(4)
    if ( phi <= 1 ):
        # lean combustion
        n[0] = 1;
        n[1] = 3.76*a_s/phi;
        n[2] = a_s/phi-1/2;
    else :    
        n[0] = 2*a_s/phi;
        n[1] = 3.76*a_s/phi;
        n[2] = 0;
        n[3] = 1-2*a_s/phi
    
    # Nombre total de moles
    N = sum(n)
    if T<1000 :
        m_r = 0;
        for i in range(len(n)):
            Y[i] = n[i]/N
            m_r = m_r + Y[i]*Mi[i];
    
        y_r = 0
        # compute residual mole fraction
        if f !=0 :
            y_r = 1/(1 + m_r/m_fa * (1/f-1));
        
        # compute total mole fractions in mixture
        for i in range(len(n)):
            Y[i] = Y[i]*y_r;
        
        # fuel mole fraction based on all moles
        y_fuel = y_fuel*(1 - y_r);
        
        # include intake N2 and O2
        Y[2] = Y[2] + y_N2*(1 - y_r);
        Y[3] = Y[3] + y_O2*(1 - y_r);
        
    if T>1000 :
        for i in range(len(n)):
            Y[i] = n[i]/N
    
    h = 0
    s = 0
    Cp = 0
    MW = 0
    if T<1000:
        # compute properties of mixture
        h = h_fuel*y_fuel;
        s = (so_fuel-np.log(max(y_fuel,1e-15)))*y_fuel;
        Cp = cp_fuel*y_fuel;
        MW = m_fuel*y_fuel;
        
    
    # compute component properties according to curve fits
    cpo = np.zeros(6);
    ho = np.zeros(6);
    so = np.zeros(6);
    for i in range(len(n)):
        cpo[i] = A[i][0] + A[i][1]*T + A[i][2]*T**2 + A[i][3]*T**3 + A[i][4]*T**4;
        ho[i] = A[i][0] + A[i][1]/2*T + A[i][2]/3*T**2 + A[i][3]/4*T**3 + A[i][4]/5*T**4 + A[i][5]/T;
        so[i] = A[i][0]*np.log(T) + A[i][1]*T + A[i][2]/2*T**2 + A[i][3]/3*T**3 +A[i][4]/4*T**4 +A[i][6];
    
    for i in range(len(n)):
        if(Y[i]>1.e-25):
            h = h + ho[i]*Y[i];
            s = s + Y[i]*(so[i]-np.log(Y[i]));
            Cp = Cp+cpo[i]*Y[i];
            MW = MW + Y[i]*Mi[i]
    
    R = 8.31434/MW;
    h = R*T*h;
    u = h-R*T;
    v = R*T/p;
    s = R*(-np.log(p/101.325)+s);
    Cp = R*Cp; 
    dvdT = v/T; 
    dvdP = -v/p; 
    return Y, h, u, s, v, R, Cp, MW, dvdT, dvdP

def Equilibrium( T, P, phi):
    """
    Computes the thermodynamic properties for the general combustion equation,
    H2 + l/2*(O2 + 3.76*N2) -> n0*H2 + n1*O2 + n2*H2O + n3*N2 + n4*O + n5*N + n6*NO + n7*OH

    Parameters
    ----------
    T : integer
        Temperature ~ [K].
    p : integer
        Pressure ~ [kPa].
    phi : integer
        Richness of the fresh gas mixture ~ [/].

    Returns
    -------
    ierr : error code :
         0 = success
         1 = singular matrix
         2 = maximal pivot error in gaussian elimination
         3 = no solution in maximum number of iterations
         4 = result failed consistency check sum(Y)=1
         5 = failure to obtain initial guess for oxygen concentration
         6 = negative oxygen concentration in initial guess calculation
         7 = maximum iterations reached in initial guess solution
         8 = temperature out of range
         9 = pressure out of range
        10 = equivalence ratio too lean
        11 = equivalence ratio too rich, solid carbon will be formed for given fuel
    Y : mole fraction of constituents
        y(0)  : H2O
        y(1)  : N2
        y(2)  : O2
        y(3)  : H2
        y(4)  : H
        y(5)  : O
        y(6)  : NO
        y(7)  : OH
    """
    Patm = 101.325
    
    # Get fuel composition information
    alpha, beta, gamma, delta, h_fuel, so_fuel, cp_fuel, m_fuel, Fs, q = fuel(T)
    
    #print(h_fuel)
    
    # initialize outputs
    Y = np.zeros(8);
    # compute cp,h,s
    # initialize h, etc to zero
    a_s = alpha + beta/4 - gamma/2;
    
    # mole fraction of fuel, O2, N2
    y_1 = 1 / (1 + 4.76*a_s/phi);  # mole fraction for one mole of reactant
    y_fuel = y_1; # assuming 1 mole fuel
    MW = 0
    Cp = 0
    h = 0
    s = 0
    u = 0;
    v = 0;
    R = 0;
    dvdT = 0;
    dvdP = 0;
    dMWdT = 0;
    dMWdP = 0;
    
    # solution parameters
    prec = 1e-3;
    MaxIter = 50 #20;
    
    # square root of pressure (used many times below)
    PATM = P/Patm;
    sqp = np.sqrt(PATM);
    
    if (T<400) or (T>4000) :
        ierr = 8;
        return ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP
    
    if (P<20) or (P>30000) :
        ierr = 9;
        return ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP
    
    if (phi<0.01):
        ierr = 10;
        return ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP
    
    
    l = 1/phi
    
    # Equilibrium constant curve fit coefficients.  
    # Valid in range: 600 K < T < 4000 K
    #         Ai           Bi            Ci            Di             Ei
    Kp = [[  0.432168,    -0.112464e+5,  0.267269e+1, -0.745744e-4,  0.242484e-8],  # 1/2 H2-H
       [  0.310805,    -0.129540e+5,  0.321779e+1, -0.738336e-4,  0.344645e-8 ],    # 1/2 O2-O
       [ -0.141784,    -0.213308e+4,  0.853461,     0.355015e-4, -0.310227e-8 ],    # 1/2 H2 + 1/2 O2-OH
       [  0.150879e-1, -0.470959e+4,  0.646096,     0.272805e-5, -0.154444e-8 ],    # 1/2 O2 + 1/2 N2 - NO
       [ -0.752364,     0.124210e+5, -0.260286e+1,  0.259556e-3, -0.162687e-7 ]]    # H2+1/2 O2-H2O
    
    K = np.zeros(len(Kp))
    for i in range(len(Kp)):
        log10ki = Kp[i][0]*np.log(T/1000) + Kp[i][1]/T  +  Kp[i][2] + Kp[i][3]*T + Kp[i][4]*T*T;
        K[i] = 10**log10ki;     # pg90 eq(3.104)
    
    c1 = K[0]/sqp
    c2 = K[1]/sqp
    c3 = K[2]
    c4 = K[3]
    c5 = K[4]*sqp
    
    ierr, y2, y3, y4 = guess( T, phi, alpha, beta, gamma, delta, c5, 0 );
    
    if ( ierr != 0 ):
        return ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP    
    
    D1 = l/beta
    D2 = 3.76*l/(2*beta)
    A = np.zeros((3,3));
    
    final = 0;
    
    for jj in range(MaxIter+1) :
        
        sqy2 = np.sqrt(y2);
        sqy3 = np.sqrt(y3);
        sqy4 = np.sqrt(y4);
    
        y5 = c1*sqy4;
        y6 = c2*sqy3;
        y8 = c3*sqy4*sqy3;
        y7 = c4*sqy2*sqy3;
        y1 = c5*sqy3*y4;
    
        d54 = 0.5*c1/sqy4;
        d63 = 0.5*c2/sqy3;
        d83 = 0.5*c3*sqy4/sqy3;
        d84 = 0.5*c3*sqy3/sqy4;
        d72 = 0.5*c4*sqy3/sqy2;
        d73 = 0.5*c4*sqy2/sqy3;
        d13 = 0.5*c5*y4/sqy3;
        d14 = c5*sqy3
    
    	# form the Jacobian matrix        
        A = [ [ 1+d72,  1+d13+d63+d73+d83,   1+d14+d54+d84 ],
          [ d72,   2+d13+d63+d83+d73-2*D1*d13-D1*d83,  d14+d84-2*D1*d14-2*D1-D1*d54-D1*d84 ],
          [ 2+d72, d73-2*D2*d13-D2*d83, -2*D2*d14-2*D2-D2*d54-D2*d84 ] ]
        
        if ( final==1 ) :
            break;
        
        B = [ -(y1+y2+y3+y4+y5+y6+y7+y8-1),
             -(y1 + 2.*y3 + y6 + y7 + y8 -2*D1*y1 - 2*D1*y4 - D1*y5 - D1*y8),
             -(2.*y2 + y7 - 2*D2*y1 - 2*D2*y4 - D2*y5 - D2*y8)]
        
        B, ierr = gauss( A, B );
        if ( ierr != 0 ) :
            return ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP
        
        y2 = y2 + B[0];
        y3 = y3 + B[1];
        y4 = y4 + B[2];
    
        nck = 0;
        if ( abs(B[0]/y2) > prec ) :
            nck = nck+1;

        if ( abs(B[1]/y3) > prec ) :
            nck = nck+1;

        if ( abs(B[2]/y4) > prec ) :
            nck = nck+1;

        if( nck == 0 ) :
            final = 1;
            continue;    
    
        if (jj>=MaxIter) :
            ierr = 3;
            return ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP
        
    Y = np.array([y1, y2, y3, y4, y5, y6, y7, y8])
    
    if( abs( sum(Y)-1 ) > 0.000001 ) :
        ierr = 4;
        return ierr, Y, h, u, s, v, R, Cp, MW, dvdT, dvdP    
    
    ierr = 0;
    return ierr, Y


def guess( T, phi, alpha, beta, gamma, delta, c5, c6 ):
    """
    First approximation of the concentration of [N2, O2 et H2]

    Parameters
    ----------
    T : integer
        Temperature ~ [K].
    phi : integer
        Equivalence ratio ~ [/].
    alpha : integer
        Number of carbon molecule.
    beta : integer
        Number of hydrogen molecules.
    gamma : integer
        Number of oxygen molecules.
    delta : integer
        Number of nitrogen molecules.
    c5 : integer
        Concentration of H2O ~ [mol_H2O/mol].
    c6 : integer
        Contration of CO2 -> set to 0.

    Returns
    -------
    ierr : error code :
         0 = success
         5 = failure to obtain initial guess for oxygen concentration
         6 = negative oxygen concentration in initial guess calculation
         7 = maximum iterations reached in initial guess solution
    y2 : integer
        Concentration of nitrogen.
    y3 : integer
        Concentration of oxygen.
    y4 : integer
        Concentration of hydrogen.

    """
    ierr = 0;
    y2 = 0;
    y3 = 0;
    y4 = 0;

    a_s = alpha + beta/4 - gamma/2;

    n = np.zeros(4)
    if ( phi <= 1 ):
        n[0] = 1;
        n[1] = 3.76*a_s/phi;
        n[2] = a_s/phi-1/2;
    else :    
        n[0] = 2*a_s/phi;
        n[1] = 3.76*a_s/phi;
        n[2] = 0;
        n[3] = 1-2*a_s/phi


    N = sum(n);

    ox = 1;
    nIterMax=40;
    for ii in range(nIterMax):
        f = 2*N*ox - gamma - (2*a_s)/phi + (alpha*(2*c6*ox**(1/2) + 1))/(c6*ox**(1/2) + 1) + (beta*c5*ox**(1/2))/(2*c5*ox**(1/2) + 2);
        if ( f < 0 ) :
            break;
        else :
            ox = ox*0.1;
            if ( ox < 1e-37 ) :
                ierr = 5;
                return ierr, y2, y3, y4

    for ii in range(nIterMax):
        f = 2*N*ox - gamma - (2*a_s)/phi + (alpha*(2*c6*ox**(1/2) + 1))/(c6*ox**(1/2) + 1) + (beta*c5*ox**(1/2))/(2*c5*ox**(1/2) + 2);
        df = 2*N - (beta*c5**2)/(2*c5*ox**(1/2) + 2)**2 + (alpha*c6)/(ox**(1/2)*(c6*ox**(1/2) + 1)) + (beta*c5)/(2*ox**(1/2)*(2*c5*ox**(1/2) + 2)) - (alpha*c6*(2*c6*ox**(1/2) + 1))/(2*ox**(1/2)*(c6*ox**(1/2) + 1)**2);
        dox = f/df;
        ox = ox - dox;
        if ( ox < 0.0 ) :
            ierr = 6;
            return ierr, y2, y3, y4

        if ( abs(dox/ox) < 0.001 ):
            break;

    if( ii == nIterMax ):
        ierr = 7;
        return ierr, y2, y3, y4

    y2 = 0.5*(delta + a_s/phi*2*3.76)/N;
    y3 = ox;
    y4 = beta/2/N/(1+c5*np.sqrt(ox));

    return ierr, y2, y3, y4

def gauss( A, B ): 
    """
    maximum pivot gaussian elimination routine adapted
    from FORTRAN in Olikara & Borman, SAE 750468,  1975
    not using built-in MATLAB routines because they issue
    lots of warnings for close to singular matrices
    that haven't seemed to cause problems in this application
    routine below does check however for true singularity

    Parameters
    ----------
    A : array
        Coefficents of the equation system.
    B : array
        Solution of the unsolved equation system.

    Returns
    -------
    B : array
        Solution of the solved equation system.
    IERQ : error code
        0 = success
        1 = singular matrix
        2 = maximal pivot error in gaussian elimination

    """
    IERQ = 0;

    for N in range(2) :                                                          
        NP1=N+1;
        BIG = abs( A[N][N]) 
        if ( BIG < 1.0e-05) :                                            
            IBIG=N;      
            for I in range (NP1,5) :
                if( abs(A[I][N]) <= BIG ):
                    continue;

                BIG = abs(A[I][N]);
                IBIG = I;

            if(BIG <= 0.) :
                IERQ=2;
                return  B, IERQ

            if( IBIG != N) :
                for J in range(N,4):                                                              
                    TERM = A[N][J];
                    A[N][J] = A[IBIG][J];
                    A[IBIG][J] = TERM                                                         

                TERM = B[N];
                B[N] = B[IBIG];
                B[IBIG] = TERM;

        for I in range(NP1,3) :  
            TERM = A[I][N]/A[N][N];

            for J in range(NP1,3):
                #print(N, I, J) 
                A[I][J] = A[I][J]-A[N][J]*TERM;

            B[I] = B[I]-B[N]*TERM;
    
    if( abs(A[2][2]) > 0.0 ):       
        B[2] = B[2]/A[2][2];
        B[1] = (B[1]-A[1][2]*B[2])/A[1][1];
        B[0] = (B[0]-A[0][1]*B[1]-A[0][2]*B[2])/A[0][0];
    else :
        IERQ=1;
        return  B, IERQ
        
    return B, IERQ