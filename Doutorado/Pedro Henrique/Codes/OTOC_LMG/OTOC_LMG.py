

import numpy as np
from numpy import linalg as LA
# from scipy.linalg import fractional_matrix_power
from scipy.linalg import expm
# from scipy.linalg import logm
# from functools import reduce
# from scipy.signal import find_peaks
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 300
dim = np.int(2*j+1)
gx = 1
vect = np.linspace(0, 20, 150)

    
    
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

# Commutator 
def Comm(A,B):
    return np.dot(A,B)-np.dot(B,A)
    

# Dagger
def dag(A):
    return np.transpose(np.conjugate(A))


#Initial LMG hamiltonian
h=0
H0 = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j))


#Initial density matrix of the system
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index] 
V0 = V0[:,index]

# beta = 1


h = 1

#Final LMG hamiltonian
Hf = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j))


############################################
#
# FOTOC => < W'(t)V'W(t)V >
#
############################################


delta = 10**(-4)

rho0 = np.dot(V0[:,0].reshape(dim,1), np.conjugate(V0[:,0]).reshape(1, dim))
# W0 = np.dot(V[:,2].reshape(dim,1), np.conjugate(V[:,2]).reshape(1, dim))
# V0 = np.dot(V[:,1].reshape(dim,1), np.conjugate(V[:,1]).reshape(1, dim))
G = np.dot(Jx(j), Jx(j))
W0 = expm(1j*delta*G)

O = np.zeros(len(vect)) #OTOC

for tt in range(len(vect)):
    
    U = expm(-1j*Hf*vect[tt])  #time evolution operator 
    
    Wt = np.dot(U, np.dot(W0, dag(U)))
    # Gt = np.dot(U, np.dot(G, dag(U)))
    
    O[tt] = np.dot(np.conjugate(V0[:,0]), np.dot(np.dot(dag(Comm(Wt, rho0)), Comm(Wt, rho0)), V0[:,0]))
    
    
plt.plot(gx*vect, O, label=r"$h_f={}, N={}$".format(h, 2*j))
plt.legend(fontsize = 12)
plt.ylabel(r"$\log(O(t))$", fontsize=17)
plt.xlabel(r"$\gamma_x t$", fontsize=17)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(-3,-3))
plt.tight_layout()



#############################################
#
# Long-time average FOTOC
#
#############################################



# delta = 10**(-4)

# V0 = np.dot(V[:,0].reshape(dim, 1), V[:,0].reshape(1, dim))
# rho0 = V0
# G = np.dot(Jx(j), Jx(j))
# W0 = expm(1j*delta*G)


# hf = np.linspace(0, 2, 50) # final magnetic fields
# T = np.linspace(0, 50, len(vect))
# O = np.zeros((50, len(vect))) #OTOC
# Oa = np.zeros(len(hf)) #time averaged OTOC


# for h in range(len(hf)):

#     #Final LMG hamiltonian
#     Hf = -hf[h]*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))
    
    
#     for tt in range(len(vect)):
        
#         U = expm(-1j*Hf*vect[tt])  #time evolution operator 
        
#         Wt = np.dot(U, np.dot(W0, dag(U)))
#         # Gt = np.dot(U, np.dot(G, dag(U)))
        
#         O[h, tt] = np.trace(np.dot(rho0, np.dot(-Comm(dag(Wt), V0), Comm(Wt, V0))))
#         # O[tt] = delta**2*(np.trace(np.dot(rho0, np.dot(Gt, Gt))) - np.trace(np.dot(rho0, Gt))**2)
        
#     Oa[h] = np.trapz(O[h,:], T)/50
    
    
# plt.plot(hf, Oa, label=r"$j=100, T=50, \delta=10^{-4}$")
# plt.legend(fontsize = 12)
# plt.ylabel(r"$\bar{O}$", fontsize=17)
# plt.xlabel(r"$h$", fontsize=17)
# # plt.ticklabel_format(axis='y', style='sci', scilimits=(-3,-3))
# plt.tight_layout()


# plt.savefig("TimeAverage_OTOC_j100.pdf")

############################################
#
# Participation ratio: < [Wt, V0]^2 >
#
############################################


# Npc = np.zeros(len(vect)) #Participation ratio


# for tt in range(len(vect)):
    
#     U = expm(-1j*vect[tt]*Hf)  #time evolution operator
    
#     Pt = 0  #sum of all squared probabilities
#     for k in range(dim):
        
#         Pt = Pt + np.abs(np.dot(V[:,k], np.dot(U, V[:,0])))**4
    
#     Npc[tt] = 1/Pt


# plt.plot(gx*vect, 1/Npc, label=r"$N_{pc}$")
# # plt.plot(gx*vect, L, label=r"$\mathcal{L}(t)$")
# # plt.xlabel(r"$\gamma_x t$", fontsize=17)
# plt.legend(fontsize=12)

# plt.savefig("Npc_LE_j100_h1.2.pdf")