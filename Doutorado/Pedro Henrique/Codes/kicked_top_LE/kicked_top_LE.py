

import numpy as np
from numpy import linalg as LA
# from scipy.linalg import fractional_matrix_power
from scipy.linalg import expm
# from scipy.linalg import logm
# from scipy.signal import find_peaks
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


p = np.pi/2
tau = 0.1
k = 0.1
j = 300
dim = np.int(2*j+1)


#Defining the Jx, Jy and Jz in z basis
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
            
            p2 = np.int(m2+j)
            p1 = np.int(m1+j)
            
            Jy[p1,p2] = 1j*(aut1-aut2)
            
            aut1 = aut2 = 0
    return Jy

def Jz(j):
    Jz = np.zeros((dim, dim))
    
    for m in np.arange(-j, j+1, 1):
        p = np.int(m+j)
        Jz[p,p] = -m
    return Jz



#Finding the eigenvalues and eigenvectors of H when k=0
E0, V0 = LA.eig((p/tau)*Jy(j))
index = np.argsort(E0) #Get the index list that would sort E0 and V0
E0 = E0[index] 
V0 = V0[:,index]


#Initializing the system in the ground state of H (k=0)
psi0 = V0[:,0]



###########################################
#
# Loshmidt echo
#
###########################################

N = 300
L = np.zeros(N)

# Floquet operator
U = np.dot(expm(-1j*(p/tau)*Jy(j)), expm(-1j*k/(2*j)*np.dot(Jz(j),Jz(j))))

for tt in range(1,N):
    
    if tt <= 2:
        
        Un = LA.matrix_power(U,tt)
        
        L[tt] = np.abs(np.dot(psi0, np.dot(Un, psi0)))
    
        U_aux = Un
        
    if tt > 2:
        Up = np.dot(U, U_aux)
        
        L[tt] = np.abs(np.dot(psi0, np.dot(Up, psi0)))

        U_aux = Up

plt.plot(np.arange(N), L, 'o-', label=r"$k={}$".format(k))
plt.legend(fontsize=10)
