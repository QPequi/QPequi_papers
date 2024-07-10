

import numpy as np
from numpy import linalg as LA
from scipy.linalg import expm
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True

j = 100
dim = np.int(2*j+1)
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



#Initial hamiltonian
alpha = 0
H0 = -4*(1-alpha)/2/j*np.dot(Jx(j), Jx(j)) + alpha*Jz(j) + alpha*j*np.identity(dim)

#Diagonalization of initial hamiltonian
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index] 
V0 = V0[:,index]


#beta = 1

#Initial density matrix
#psi0 = expm(-beta*H0)/np.trace(expm(-beta*H0))
# psi0 = np.dot(V0[:,0].reshape(dim,1), V0[:,0].reshape(1,dim))

alpha = 1.2
#Final hamiltonian

Hf = -4*(1-alpha)/2/j*np.dot(Jx(j), Jx(j)) + alpha*Jz(j) + alpha*j*np.identity(dim)

# Diagonalization of the final hamiltonian
Ef, Vf = LA.eig(Hf)
index = np.argsort(Ef) #Get the index list that would sort E and V
Ef = Ef[index] 
Vf = Vf[:,index]


########################################################
#
# Diagonal entropy  S = -sum rho_nn*log(rho_nn)
#
########################################################


DE = np.zeros(len(vect)) #Diagonal entropy

for tt in range(len(vect)):
    U = expm(-1j*vect[tt]*Hf)
    
    for n in range(dim):
        DE[tt] = DE[tt] - np.abs(np.dot(V0[:,n], np.dot(U, V0[:,0])))**2*np.log(np.abs(np.dot(V0[:,n], np.dot(U, V0[:,0])))**2)
    

    
plt.plot(vect, DE, label="$h={}, j={}$".format(alpha, j))
plt.legend(fontsize = 12)
# plt.xlabel("$\gamma_x t$", fontsize=17)
# plt.xlim([0, 3])
plt.ylabel("$S_d(t)$", fontsize=17)
plt.tight_layout()
   

#plt.savefig("DiagEntropy_LMG_1.pdf")



    
    
    
    
    