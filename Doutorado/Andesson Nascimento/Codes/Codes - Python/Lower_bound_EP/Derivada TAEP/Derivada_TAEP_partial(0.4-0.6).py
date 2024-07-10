import os
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
import math as mt
import time

initial_t = time.time()
# Calculo da derivada da média temporal de produção de entropia

j = 500
dim = int(2*j+1)
g = 0.5
vect = np.linspace(0,30,200)
hvec = np.linspace(0.4,0.6,50) # Be carreful! h =/= 0 because LL becomes indeterminate

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

L = np.zeros(len(vect))  # Loschmidt Echo
LL = np.zeros(len(vect)) # Bures Angle
EP = np.zeros(len(vect)) # Minimal Entropy Production - 1 term
TAEP = np.zeros(len(hvec)) #  Time Average - Entropy production
dTAEP = np.zeros(len(vect))

for h in range(len(hvec)):
#Final hamiltonian
    H = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*hvec[h]*Jz(j)
    for tt in range(len(vect)):
        U = expm(-1j*H*vect[tt])
        L[tt] = round(np.abs(np.dot(dag(V0[:,0]), np.dot(U, V0[:,0])))**2,10) 
        #if L[tt] <= 0:
        #    L[tt] = 1e-16
        LL[tt] = mt.acos(L[tt])
        EP[tt] = 8/(mt.pi**2)*LL[tt]**2
    TAEP[h] = (1/np.max(vect))*np.trapz(EP, vect)
    
dTAEP = np.gradient(TAEP)
dh = np.gradient(hvec)   

Der = dTAEP/dh


# # Salving data in .txt 
# np.savetxt(os.path.join(os.getcwd(),'Der_TAEP_j={}.txt'.format(j), Der))

# Plotting Time Average Entropy Production

plt.plot(hvec,Der)
plt.ylabel(r'$\frac{d\bar{\langle \Sigma \rangle}}{dh}$', fontsize=14)
plt.xlabel('h', fontsize=14)
plt.axvline(0.5, color="red",linestyle='--',label='Dynamical Critical Point')
plt.legend()

fig = plt.gcf()
fig.set_size_inches(8, 6)
fig.set_dpi(300)

# Saving figure in .SVG
plt.savefig('DerTAEP_part_j={}.svg'.format(j), format='svg')

final_t  = time.time()

total_t = (final_t - initial_t)*0.0167

print(f"Tempo de Execução = {total_t} Min")










