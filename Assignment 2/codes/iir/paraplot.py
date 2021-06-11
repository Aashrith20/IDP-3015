# PLOTS OF THE LOWPASS CHEBYSCHEV FILTER OF ORDER N AND
# 0.3184 < epsilon < 0.6197
import numpy as np
import matplotlib.pyplot as plt

plt.figure()
N = 4
epsVals = np.arange(0.35,0.61,0.05)
for epsilon in epsVals:
    Omega1 = np.arange(0,1,0.01)
    H1 = 1/np.sqrt(1 + epsilon**2*np.cos(N*np.arccos(Omega1))**2)
    Omega2 = np.arange(1,3.01,0.01)
    H2 = 1./np.sqrt(1 + (epsilon**2)*np.power(np.cosh(N*np.arccosh(Omega2)),2))
    H = np.concatenate((H1,H2))
    Omega = np.concatenate((Omega1,Omega2))
    plt.plot(Omega,H,label='$\epsilon$ = '+str(round(epsilon,2)))

plt.xlabel('$\Omega$')
plt.ylabel('$|H_{a,LP}(j\Omega)|$')
plt.legend()
plt.grid()
plt.show()

