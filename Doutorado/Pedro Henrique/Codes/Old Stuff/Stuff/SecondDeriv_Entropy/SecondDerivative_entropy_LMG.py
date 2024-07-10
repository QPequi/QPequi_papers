

import numpy as np
from numpy import linalg as LA
from scipy.linalg import expm
from scipy.linalg import logm
import matplotlib.pyplot as plt

j=100
dim = np.int(2*j+1)
gx=1
gy = 0
vect = np.linspace(0,30,300)


#Polar Coordinates
theta = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)

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
    

#Initial hamiltonian
h=0
H0 = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))


#Initial density matrix of the system
E, V = LA.eig(H0)
index = np.argsort(E) #Get the index list that would sort E and V
E = E[index] 
V = V[index]

psi0 = np.dot(V[:,0].reshape(dim,1), np.conjugate(V[:,0]).reshape(1,dim))  #density matrix

#Final hamiltonian
h = 2
Hf = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))


##############################################
#
# Second derivative of the entropy
#
##############################################


SDEntro = np.zeros(len(vect)) #Second derivative of the entropy

for tt in range(len(vect)):
    U = expm(-1j*Hf*vect[tt])
    Ut = np.transpose(np.conjugate(U))
    SDEntro[tt] = np.trace(np.dot(psi0, Comm(H0, np.dot(np.dot(U, Comm(H0, logm(psi0))), Ut))))
    
    
plt.plot(gx*vect, np.real(SDEntro), label="$h={}, \gamma_x={}, j={}$".format(h,gx,j))
plt.xlabel("$\gamma_x t$", fontsize=17)
plt.ylabel(r"$J(\rho (t),t)$", fontsize=17)
plt.legend()
plt.tight_layout()
#plt.savefig("SecondDeriv_Entropy_j300.pdf")





