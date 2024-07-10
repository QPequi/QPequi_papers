

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 500
dim = np.int(2*j+1)
vect = np.linspace(0, 100, 150)



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
alpha = 0.4
H0 = -2*(1-alpha)/j*np.dot(Jx(j), Jx(j)) + alpha*(Jz(j) + j*np.identity(dim))


#Initial hamiltonian's eigenvalues and eigenvectors
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index] 
V0 = V0[:,index]


#initial state
psi0 = V0[:,0]  #pure ground-state of H0

n = 25
Npc_mean = np.zeros(n)
i = 0

Npc = np.zeros(len(vect)) #Participation ratio
beta = np.linspace(0, 3, n)

for b in beta:
    
    Hf = H0 + b*(Jz(j) + j*np.identity(dim))
    
    Npc = np.zeros(len(vect))
    
    for tt in range(len(vect)):
        
        U = expm(-1j*vect[tt]*Hf)  #time evolution operator
        
        Pt = 0  #sum of all squared probabilities
        for k in range(dim):
            
            Pt = Pt + np.abs(np.dot(np.conjugate(V0[:,k]), np.dot(U, V0[:,0])))**4
        
        Npc[tt] = 1/Pt
    
    Npc_mean[i] = (1/max(vect))*np.trapz(Npc, vect)
    i = i+1

plt.plot(beta, Npc_mean, label=r"$N={}$".format(2*j))
plt.axvline(1, c='red', linestyle='--')
plt.xlabel(r"$\lambda$", fontsize=18)
plt.ylabel(r"$\overline{N_{pc}}$", fontsize=18)
plt.legend(fontsize=15)
plt.tight_layout()
# plt.savefig("Participation_N1000_ParametLMG.pdf")
    




# plt.plot(gx*vect, tau, label=r"$\tau^{R}_{\alpha}(\rho(t)||\rho_0)$")
# plt.plot(gx*vect, 20*L, label=r"$\mathcal{L}(t)$")
# plt.xlabel(r"$\gamma_x t$", fontsize=17)
# plt.title(r"$h={},\alpha={},j={}$".format(h,alfa,j))
# plt.legend(fontsize=12)
# plt.tight_layout()
# plt.savefig("QSL_LE_Renyi_purestate_j300.pdf")

