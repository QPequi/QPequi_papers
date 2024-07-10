

import numpy as np
from scipy.linalg import expm
from scipy import linalg as LA
import matplotlib.pyplot as plt
from scipy.signal import find_peaks #it helps to identify t_c (peaks)


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True



j = 100
dim = int(2*j+1)
# omega = 0


#Defining the Jx, Jy and Jz in z basis
def Jx(j):
    Jx = np.zeros((dim, dim))
    
    aut1=0
    aut2=0
    
    for m1 in np.arange(-j, j+1, 1): #row
        for m2 in np.arange(-j, j+1, 1): #column
            if m1 == m2+1:
                aut1 = np.sqrt(j*(j+1)-m2*(m2+1))/2
            if m1 == m2-1:
                aut2 = np.sqrt(j*(j+1)-m2*(m2-1))/2
            
            p1 = int(m1+j)
            p2 = int(m2+j)
            
            Jx[p1,p2] = aut1+aut2
            
            aut1 = aut2 = 0
    return Jx

def Jy(j):
    Jy = np.zeros((dim, dim),dtype=complex)
    
    aut1=0
    aut2=0
    
    for m1 in np.arange(-j, j+1, 1): #row
        for m2 in np.arange(-j, j+1, 1): #column
            if m1 == m2+1:
                aut1 = np.sqrt(j*(j+1)-m2*(m2+1))/2
            if m1 == m2-1:
                aut2 = np.sqrt(j*(j+1)-m2*(m2-1))/2
            
            p2 = int(m2+j)
            p1 = int(m1+j)
            
            Jy[p1,p2] = 1j*(aut1-aut2)
            
            aut1 = aut2 = 0
    return Jy

def Jz(j):
    Jz = np.zeros((dim, dim))
    
    for m in np.arange(-j, j+1, 1):
        p = int(m+j)
        Jz[p,p] = -m
    return Jz

def H_LMG(h):
    return -(1/(2*j))*np.dot(Jz(j),Jz(j)) - h*Jx(j)


Energies = np.zeros((dim, 20))

i=0
for h in np.linspace(0,1,20):
    
    
    #Pre-quench eigenvalues and eigenvectors
    E, V = LA.eig(H_LMG(h))
    index = np.argsort(E) #Get the index list that sorts E and V
    Energies[:,i] = E[index]
    
    i += 1
    # plt.subplot(2,2,i)
    # plt.hist(E, bins=80, density=True, color='#003366', label=r"$\Omega={}$".format(omega))
    # plt.legend(fontsize=10, frameon=False)
    # plt.ylim([0,10])
    
    # if i==1 or i==3:
    #     plt.ylabel(r"$\nu(E)$")
    
    # if i==3 or i==4:
    #     plt.xlabel(r"$E_n/N$")
    
    # i += 1
    

for i in range(5):
    plt.plot(np.linspace(0,1,20), Energies[i,:])
    
plt.tight_layout()
# plt.savefig("DOS_LMG_N2000_2.pdf")










