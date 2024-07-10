
import numpy as np
from numpy import linalg as LA
from scipy.linalg import expm
#from scipy import integrate
import matplotlib.pyplot as plt

j=1
dim = np.int(2*j+1)
gx=1
gy = 0
vect = np.linspace(0,100,300)


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


#Initial hamiltonian
h=0
H0 = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))


#Diagonalization of the inicial hamiltonian
E, V = LA.eig(H0)
index = np.argsort(E) #Get the index list that would sort E and V
E = E[index] 
V = V[index]

#density matrix
psi0 = np.dot(V[:,0].reshape(dim,1), np.conjugate(V[:,0]).reshape(1,dim))  


#Matrix of spin coherent state
Omega = np.zeros((len(theta),len(phi)), dtype=object)
state_spin = V[:,dim-1].reshape(dim,1)  #Get the eigenvector associated with the max eigv

for gg in range(len(theta)):
    for pp in range(len(phi)):
        #Generate the spin coherent states
        Omega[gg,pp] = np.dot(expm(-1j*theta[gg]*Jz(j)), np.dot(expm(-1j*phi[pp]*Jy(j)),state_spin))
        
        
#Final hamiltonian
h = 0.8
Hf = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))


#########################################################
#
# ENTROPY
#
#########################################################

Entro = np.zeros(len(vect))

#Evolved state
for tt in range(len(vect)):
    U = expm(-1j*vect[tt]*Hf)
    Ut = np.transpose(np.conjugate(U))
    psit = np.dot(U, np.dot(psi0, Ut))  #density matrix evolved
    
    sq = np.zeros((len(theta),len(phi))) #Stores sin(theta)*Q*log(Q) to perform the integral later
    
    E_theta = np.zeros(len(theta)) #Stores the integration just over phi
    for gg in range(len(theta)):
        for pp in range(len(phi)):
            Q = np.real(np.dot(np.transpose(np.conjugate(Omega[gg,pp])), np.dot(psit, Omega[gg,pp])))
            
            if Q == 0: #Set 0*log(0) = 0
                sq[gg,pp] = 0
            else: #Sin(theta) comes from the jacobian 
                sq[gg,pp] = -np.sin(theta[gg])*Q*np.log(Q) 
        
        E_theta[gg] = np.trapz(sq[gg,:], phi)
    Entro[tt] = dim*np.trapz(E_theta, theta)/4/np.pi


plt.plot(vect, Entro)

















        