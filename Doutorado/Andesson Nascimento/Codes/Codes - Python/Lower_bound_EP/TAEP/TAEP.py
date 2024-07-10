import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
import math as mt
import os

# This code computes the Time Average Entropy Production (TAEP) to LMG model with hamiltonian 
#   H0 = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h0*Jz(j) with J = 1 (and dynamical critical point h = 1)
# To compute TAEP, I did a variation on parameter h from 0 to 1,5, crossing the dynamic critical point.

vect = np.linspace(0,30,150)
hvec = np.linspace(0.1,0.8,50) # Be carreful! h =/= 0 because LL becomes indeterminate
lvec = [1] # System's Dimension

L = np.zeros(len(vect))  # Loschmidt Echo
LL = np.zeros(len(vect)) # Bures Angle
EP = np.zeros(len(vect)) # Minimal Entropy Production - 1 term
s_ub = np.zeros(len(vect)) # Upper Bond of Entropy production
TAEP = np.zeros(len(hvec)) # Time Average - Minimal Entropy production
TAEP_s = np.zeros(len(hvec)) # Time Average - Upper bound Entropy production

for ll in range(len(lvec)):
    j = lvec[ll]*100
    dim = int(2*j+1)
    g = 0.5

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

#####################################################################
# Initial conditions
#####################################################################
#Initial hamiltonian
    h0 = 0
    H0 = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h0*Jz(j) 
    #Initial density matrix of the system
    E0, V0 = LA.eig(H0)
    index = np.argsort(E0) #Get the index list that would sort E and V
    E0 = E0[index]
    V0 = V0[:,index]

#####################################################################
# TAEP
#####################################################################

    for h in range(len(hvec)):
        #Final hamiltonian
        H = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*hvec[h]*Jz(j)
        for tt in range(len(vect)):
            U = expm(-1j*H*vect[tt])
            L[tt] = round(np.abs(np.dot(dag(V0[:,0]), np.dot(U, V0[:,0])))**2,12) 
            if L[tt] <= 0:
                L[tt] = 1e-16
            LL[tt] = mt.acos(L[tt])
            s_ub[tt] = - np.log(1-(2/mt.pi)*LL[tt])
            EP[tt] = (8/mt.pi**2)*LL[tt]**2
        TAEP[h] = (1/np.max(vect))*np.trapz(EP, vect)
        TAEP_s[h] = (1/np.max(vect))*np.trapz(s_ub, vect)
            # Salving data in .txt 
            # np.savetxt(os.path.join(os.getcwd(),'TAEP_j={}.txt'.format(j)), TAEP)
            
# Plotting Time Average Entropy Production

    # plt.plot(hvec,TAEP,label='j={}'.format(lvec[ll]*100))
    # plt.ylabel(r'$\bar{\langle \Sigma \rangle}$', fontsize=14)
    # plt.xlabel('h', fontsize=14)
    # plt.axvline(0.5, color="red",linestyle='--')
    # plt.legend()


