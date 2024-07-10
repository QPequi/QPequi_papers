


import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm



plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True

n = 6
SDev = np.zeros(n)
Jvalues = np.arange(100, 100*(n+1), 100)


for hf in [0.4, 0.6, 0.8, 1, 1.2]:
    
    i = 0
    for j in Jvalues:

        J = 1
        g = 0.5
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
        
        # Commutator 
        def Comm(A,B):
            return np.dot(A,B)-np.dot(B,A)
        
        
        
        #Initial LMG hamiltonian
        h=0
        H0 = -h*Jz(j) - (J/2/j)*np.dot(Jx(j),Jx(j)) - (J*g/2/j)*np.dot(Jy(j),Jy(j))
        
        
        #Initial hamiltonian's eigenvalues and eigenvectors
        E0, V0 = LA.eig(H0)
        index = np.argsort(E0) #Get the index list that would sort E and V
        E0 = E0[index] 
        V0 = V0[:,index]
        
        
        ###############################################
        #
        # Stantard deviation of the LMG model
        # 
        ###############################################
        
        # Quenched LMG homailtonian
        
        H = -hf*Jz(j) - (J/2/j)*np.dot(Jx(j),Jx(j)) - (J*g/2/j)*np.dot(Jy(j),Jy(j))
        
        SDev[i] = np.sqrt(np.dot(np.conjugate(V0[:,0]), np.dot(np.dot(H,H), V0[:,0])) - np.dot(np.conjugate(V0[:,0]), np.dot(H, V0[:,0]))**2)
        
        i = i+1
        
        
    plt.plot(2*Jvalues, SDev, 'o-', label=r"$h_f={}$".format(hf))
    plt.legend(fontsize=10)
    plt.title(r"$J = 1, \gamma = 0.5$")
    plt.xlabel(r"$N$", fontsize=15)
    plt.ylabel(r"$\Delta H_f$", fontsize=15)
    plt.tight_layout()


# plt.savefig("StandDeviation_N_AnisotropicLMG1.pdf")
    
















 