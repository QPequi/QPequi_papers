
import numpy as np
from scipy.linalg import expm
from scipy import linalg as LA
import matplotlib.pyplot as plt
# from scipy.signal import find_peaks #it helps to identify t_c (peaks)
from scipy.linalg import fractional_matrix_power


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 50
dim = np.int(2*j+1)
gx = 1
vect = np.linspace(0, 30, 200)



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


#Commutator [A,B]
def Comm(A,B):
    return np.dot(A,B)-np.dot(B,A)


def dag(A):
    return np.transpose(np.conjugate(A))

# Energies of LMG hamiltonian
def E(J, mz, h):
    return -2*h*mz - (gx/J)*(J*(J+1) - mz**2) + gx


#Initial LMG hamiltonian
h = 0
H0 = -2*h*Jz(j) - (gx/j)*(np.dot(Jx(j),Jx(j)) + np.dot(Jy(j),Jy(j))) + gx


beta = 2

# Z0 = 0
# for m in range(-j, j+1, 1):
#     Z0 = Z0 + np.exp(-beta*E(j,m,h))


# Diagonalizing Jz to construct Dicke states
Ez, Vz = LA.eig(Jz(j))
index = np.argsort(Ez) #Get the index list that would sort E and V
Ez = Ez[index]
Vz = Vz[:,index]


# Initial density matrix: ground state/thermal state
# rho0 = np.dot(V0[:,0].reshape(dim,1), V0[:,0].reshape(1,dim))
# rho0 = np.zeros((dim, dim))
# for n in range(-j, j+1, 1):
#     rho0 = rho0 + (1/Z0)*np.exp(-beta*E(j, n, h))*np.dot(Vz[:,n+j].reshape(dim,1), Vz[:,n+j].reshape(1,dim))


rho0 = expm(-beta*H0)/np.trace(expm(-beta*H0))

h = 0.8

#Final LMG hamiltonian
Hf = -2*h*Jz(j) - (gx/j)*(np.dot(Jx(j),Jx(j)) + np.dot(Jy(j),Jy(j))) + gx

# Zf = 0
# for m in range(-j, j+1, 1):
#     Zf = Zf + np.exp(-beta*E(j,m,gx))


# rhof_eq = np.zeros((dim, dim))
# for n in range(-j, j+1, 1):
#     rhof_eq = rhof_eq + (1/Zf)*np.exp(-beta*E(j, n, gx))*np.dot(Vz[:,n+j].reshape(dim,1), Vz[:,n+j].reshape(1,dim))
# rhof_eq = expm(-beta*Hf)/np.trace(expm(-beta*Hf))


# Final hamiltonian's eigenvalues and eigenvectors
Ef, Vf = LA.eig(Hf)
index = np.argsort(Ef) #Get the index list that would sort E and V
Ef = Ef[index]
Vf = Vf[:,index]


alpha = 0.8

RRE = np.zeros(len(vect))

for tt in range(len(vect)):
    for k in range(-j, j+1, 1):
        for n in range(-j, j+1, 1):
            for m in range(-j, j+1, 1):
                rhot = np.exp(-1j*E(j, n, h)*vect[tt])*np.exp(+1j*E(j, m, h)*vect[tt])*np.exp(-beta*E(j, k, 0))*np.dot(Vf[:,n+j], Vz[:,k+j])*np.dot(Vz[:,k+j], Vf[:,m+j])*np.dot(Vf[:,n+j].reshape(dim,1), Vf[:,m+j].reshape(1,dim))
    
    RRE[tt] = 1/(alpha-1)*np.log(np.trace(np.dot(fractional_matrix_power(rhot, alpha), fractional_matrix_power(rho0, 1-alpha))))
    
    
plt.plot(vect, RRE)
# plt.plot(vect, LE_a)
    
    
    
    
    
    
    
    
    
    
    