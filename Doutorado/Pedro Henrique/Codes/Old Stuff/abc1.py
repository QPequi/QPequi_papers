


import numpy as np
from scipy.linalg import expm
from scipy import linalg as LA
import matplotlib.pyplot as plt
# import latexify




plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 50
dim = int(2*j+1)



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


def H_LMG(chi, omega, w):
    return -(chi/(2*j))*np.dot(Jz(j),Jz(j)) - omega*Jx(j) - w*Jz(j)

H_LMG

n = dim


Energies = np.zeros((dim, n))

i = 0
for omega in np.linspace(0, 0.3, n):

    chi=1
    w = 0
    
    #Pre-quench eigenvalues and eigenvectors
    E, V = LA.eig(H_LMG(chi, omega, w))
    index = np.argsort(E) #Get the index list that sorts E and V
    E = np.real(E[index])
    V = V[:,index]

    Energies[:,i] = E
    i += 1

for j in range(5):
    
    plt.plot(np.linspace(0, 1, n), Energies[:,j])

    
    
    
    
    
    
    
