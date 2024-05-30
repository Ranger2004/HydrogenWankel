import matplotlib.pyplot as plt

from TwoZoneCycle import homogeneous

#%% PARAMETERS TO SET
phi = 1
p1 = 101.325
T1 = 350
RPM = 3000
BLOWBY = 0
thetas = 0
thetad = 70
f = 0.1

#%% RESULTS PLOT
VARIABLES = [phi, p1, T1, RPM, BLOWBY, thetas, thetad, f]
theta, vol, p, Tb, Tu, x, IMEP, NOx = homogeneous(VARIABLES)
for i in range(len(theta)):
    if Tb[i]==0:
        Tb[i]=None
    if Tu[i]==0:
        Tu[i]=None


print("""Results from simulation :
1. IMEP = {0} [kPa]
2. NOx = {1} [ppm]""".format(IMEP,NOx))
fig, axs = plt.subplots(2,2,figsize=(15, 7))

axs[0,0].plot(theta, p*10**(-2))
axs[1,0].plot(theta, Tu)
axs[1,0].plot(theta, Tb)
axs[0,1].plot(theta, vol*10**6)
axs[1,1].plot(theta, x)

axs[0,0].set(xlabel='$\\theta$ [degree]', ylabel='Pressure [bar]')
axs[1,0].set(xlabel='$\\theta$ [degree]', ylabel='Temperature [K]')
axs[0,1].set(xlabel='$\\theta$ [degree]', ylabel='Volume [$cm^3$]')
axs[1,1].set(xlabel='$\\theta$ [degree]', ylabel='Burned rate fraction')
