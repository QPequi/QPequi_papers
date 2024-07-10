import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
import math as mt

ll = 5
j = ll*100
dim = int(2*j+1)
g = 0.5
vect = np.linspace(0,30,300)
hvec = np.linspace(0.1,1,4) # post-quench parameter variation

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
def dag(A):
    return np.transpose(np.conjugate(A))

#Initial hamiltonian

h0 = 0
H0 = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h0*Jz(j)

#Initial density matrix of the system
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index]
V0 = V0[:,index]

#####################################################################
# LL and EP
#####################################################################

L = np.zeros(len(vect))  # Loschmidt Echo
LL = np.zeros(len(vect)) # Bures Angle
EP = np.zeros(len(vect)) # Minimal Entropy Production - 1 term

# for h in range(len(hvec)):
#Final hamiltonian
h = 1
H = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h*Jz(j)
for tt in range(len(vect)):
    U = expm(-1j*H*vect[tt])
    L[tt] = round(np.abs(np.dot(dag(V0[:,0]), np.dot(U, V0[:,0])))**2,10)
    LL[tt] = np.arccos(L[tt])
    EP[tt] = 8/(mt.pi**2)*LL[tt]**2

plt.plot(vect,EP,label='h=1') # Plotting # choose the graphic
#plt.plot(vect,LL, label='h={}'.format(hvec[h])) # Plotting 
plt.ylabel(r'$\langle \Sigma \rangle (t)$', fontsize=14)
plt.xlabel('t', fontsize=14)
plt.legend()

