

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True



N = 50
h0 = 0
hf = 0.01
beta = 10

vect = np.linspace(0, 2, 200)

# pseudomomenta set
K = np.array([np.pi*(2*n-1)/N for n in range(1, np.int(N/2)+1, 1)])


# Bogoliubov angle
def phi(k, g):
    return np.arctan(np.sin(k)/(g - np.cos(k)))

# Energies
def En(k, g):
    return 2*np.sqrt(np.sin(k)**2 + (g - np.cos(k))**2)


def C(k,s,t):  # s: signal + or -
    return np.cos((phi(k,hf)-phi(k,h0))/2)**2*np.exp(s*1j*t*En(k,hf))


def S(k,s,t):  # s: signal + or -
    return np.sin((phi(k,hf)-phi(k,h0))/2)**2*np.exp(s*1j*t*En(k,hf))



Z0 = 2**N

for k in K:
    Z0 = Z0*np.cosh(beta*En(k, h0)/2)**2
    


#####################################################
#
# Characteristic function Chi(t)
#
#####################################################


Chi = (1/Z0)*np.ones(len(vect))

for tt in range(len(vect)):
    
    for k in K:
        A = np.exp((1j*vect[tt] + beta)*En(k, h0))*(C(k,-1,vect[tt]) + S(k,1,vect[tt]))
        B = np.exp(-(1j*vect[tt] + beta)*En(k, h0))*(C(k,1,vect[tt]) + S(k,-1,vect[tt]))
        Chi[tt] = Chi[tt]*(A+B+2)
    

# Rate function
R = -(1/N)*np.log(np.abs(Chi)**2)
# tc, _ = find_peaks(R, height=0)

plt.plot(vect, R, label=r"$\beta={}$".format(beta))

plt.legend(fontsize=10)
plt.xlabel(r"$t$", fontsize=17)
plt.ylabel(r"$r(t)$", fontsize=17)
plt.tight_layout()
    
    
# plt.savefig("RateFunction_Ising.pdf")
    

    
    
    
    
    
