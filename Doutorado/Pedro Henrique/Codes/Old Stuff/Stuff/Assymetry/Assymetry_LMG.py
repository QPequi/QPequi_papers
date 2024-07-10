

import numpy as np
from numpy import linalg as LA
from scipy.linalg import fractional_matrix_power
from scipy.linalg import expm
from scipy.linalg import logm
from scipy.signal import find_peaks
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 100
dim = np.int(2*j+1)
gx = 1
gy = 0
vect = np.linspace(0, 100, 200)



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
            
            p1 = np.int(m1+j)
            p2 = np.int(m2+j)
            
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

#Commutator 
def Comm(A,B):
    return np.dot(A,B)-np.dot(B,A)

#Dagger
def dag(A):
    return np.transpose(np.conjugate(A))

# Trace norm
def tnorm(A):
    return np.trace(fractional_matrix_power(np.dot(dag(A), A), 1/2))
    

#Initial LMG hamiltonian
h=0
H0 = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))


#Initial density matrix of the system
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index] 
V0 = V0[:,index]



#Initial density matrix: ground state/thermal state
rho0 = np.dot(V0[:,0].reshape(dim, 1), V0[:,0].reshape(1, dim))



h = 1.2

#Final LMG hamiltonian
Hf = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))


##############################################
#
# Assymetry measure:  FL = || [\rho, L] ||_1
#
##############################################


FL = np.zeros(len(vect))  #assymetry measure

for tt in range(len(vect)):
    
    U = expm(-1j*vect[tt]*Hf)
    
    rhot = np.dot(U, np.dot(rho0, dag(U)))
    
    FL[tt] = tnorm(Comm(rhot, Jy(j)))



plt.plot(vect, FL, label=r"$h = {}, j={}$".format(h,j))
plt.xlabel(r"$t$", fontsize=17)
plt.ylabel(r"$F_{J_z}(\rho(t))$")
plt.legend(fontsize=10)
plt.tight_layout()

# plt.savefig("Assymetry_Jz_LMG.pdf")








