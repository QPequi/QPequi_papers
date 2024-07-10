

import numpy as np
from numpy import linalg as LA
from scipy.linalg import expm
import matplotlib.pyplot as plt


# plt.rcParams["font.family"] = "stix"
# plt.rcParams["mathtext.fontset"] = "cm"
# plt.rcParams["font.size"] = 15


j = 100
dim = np.int(2*j+1)
gx=1
gy = 0
vect = np.linspace(0,30,300)



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


# #Initial density matrix of the system
# E, V = LA.eig(H0)
# index = np.argsort(E) #Get the index list that would sort E and V
# E = E[index] 
# V = V[:,index]


beta = 10

#Initial density matrix
psi0 = expm(-beta*H0)/np.trace(expm(-beta*H0))


#Final hamiltonian
h = 0.8
Hf = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))



######################################################
#
# Fluctuations velocity
#
######################################################

#We calculate ||[i*rho(t), H]||_1 instead of ||d rho(t)/dt||_1

comm_norm = np.zeros(len(vect))

for tt in range(len(vect)):
    U = expm(-1j*vect[tt]*Hf)
    Ut = np.transpose(np.conjugate(U))
    
    psit = np.dot(U, np.dot(psi0,Ut))  #time evolution of psi0
    E_comm, V_comm = LA.eig(1j*Comm(psit, Hf))  #diagonalization of [rho(t), H]
    
    comm_norm[tt] = np.sum(np.abs(E_comm))  #1-norm of i*[rho(t), H]


vs = comm_norm/2

plt.plot(vect, vs, label="$h={}, \gamma_x={}, j={}$".format(h,gx,j))
plt.xlabel("$\gamma_x t$", fontsize=17)
plt.ylabel("$v_s(t)$", fontsize=17)
plt.legend()
plt.gcf().subplots_adjust(left=0.2, bottom=0.2)
#plt.savefig("FluctVelocity_LMG_j300_Hf_h14.pdf")









