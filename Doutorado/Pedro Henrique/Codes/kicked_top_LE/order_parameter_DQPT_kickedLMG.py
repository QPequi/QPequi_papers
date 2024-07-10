


import numpy as np
from numpy import linalg as LA
# from scipy.linalg import fractional_matrix_power
from scipy.linalg import expm
# from scipy.linalg import logm
# from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy import integrate


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 100
dim = np.int(2*j+1)



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

def dag(A):
    return np.transpose(np.conjugate(A))

#Initial hamiltonian
h = 1
H = -(1-h)*Jz(j) - h/2/j*np.dot(Jx(j),Jx(j))

#Initial hamiltonian's eigenvalues and eigenvectors
E0, V0 = LA.eig(H)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index] 
V0 = V0[:,index]


psi0 = V0[:,0]



######################################
#
# DQPT
#
######################################


vect = np.linspace(0, 200, 150)
h_values = np.linspace(0, 1, 100)

M = np.zeros(len(vect))
M_mean = np.zeros(len(h_values))

i = 0
for h in h_values:
    H = -(1-h)*Jz(j) - h/2/j*np.dot(Jx(j),Jx(j))
    
    for tt in range(len(vect)):
        
        U = expm(-1j*vect[tt]*H)
        psit = np.dot(U, psi0)
        
        M[tt] = np.dot(dag(psit), np.dot(Jx(j), psit))/j
        
    
    M_mean[i] = (1/np.max(vect))*integrate.simps(M, vect)
    i = i+1
    
# plt.plot(vect, Mx)
plt.plot(h_values, M_mean, label=r"DQPT $j={}$".format(j))
plt.legend(fontsize=12)
plt.xlabel(r"$s$", fontsize=17)
plt.ylabel(r"$\overline{\langle J_x\rangle/j}$", fontsize=17)
# plt.grid(True)
plt.tight_layout()

# plt.savefig("MeanMagnet_Jx_vs_h_j100.pdf")













