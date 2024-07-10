



import numpy as np
from scipy.linalg import fractional_matrix_power
from scipy.linalg import expm
from scipy.linalg import logm
# from scipy.signal import find_peaks
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 50
dim = np.int(2*j+1)
gx = 1
beta = 1
vect = np.linspace(0, 100, 300)



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
    

#Schatten 2-norm
def norm2(A):
    return np.sqrt(np.trace(np.dot(np.transpose(np.conjugate(A)), A)))


# Eigenenergies 
def E(J, jz, h):
    return -2*h*jz - (gx/J)*(J*(J+1) - jz**2) + gx



#Initial LMG hamiltonian
h = 0
H0 = -2*h*Jz(j) - (gx/j)*(np.dot(Jx(j),Jx(j)) + np.dot(Jy(j),Jy(j))) + gx*np.eye(dim)


# #Initial density matrix of the system
# E, V = LA.eig(H0)
# index = np.argsort(E) #Get the index list that would sort E and V
# E = E[index] 
# V = V[:,index]



Z0 = 0
for m in range(-j, j+1, 1):
    Z0 = Z0 + np.exp(-beta*E(j,m,h))


#Diagonalizing Jz to construct Dicke states
Ez, Vz = LA.eig(Jz(j))
index = np.argsort(Ez) #Get the index list that would sort E and V
Ez = Ez[index]
Vz = Vz[:,index]


#Initial density matrix: ground state/thermal state
# psi0 = np.dot(V[:,0].reshape(dim,1), V[:,0].reshape(1,dim))
rho0 = np.zeros((dim, dim))
for n in range(-j, j+1, 1):
    rho0 = rho0 + (1/Z0)*np.exp(-beta*E(j, n, h))*np.dot(Vz[:,n+j].reshape(dim,1), Vz[:,n+j].reshape(1,dim))


E0, V0 = LA.eig(rho0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index] 
V0 = V0[:,index]


h = 0.8

#Final LMG hamiltonian
Hf = -2*h*Jz(j) - (gx/j)*(np.dot(Jx(j),Jx(j)) + np.dot(Jy(j),Jy(j))) + gx*np.eye(dim)


Zf = 0
for m in range(-j, j+1, 1):
    Zf = Zf + np.exp(-beta*E(j,m,h))


# #  Final thermal state
# rhof_eq = np.zeros((dim, dim))
# for n in range(-j, j+1, 1):
#     rhof_eq = rhof_eq + (1/Zf)*np.exp(-beta*E(j, n, h))*np.dot(Vz[:,n+j].reshape(dim,1), Vz[:,n+j].reshape(1,dim))

rhof_eq = expm(-beta*Hf)/np.trace(expm(-beta*Hf))

####################################################
#
# Relative entropies
#
####################################################


# RE = np.zeros(len(vect)) #usual relative entropy
# RRE = np.zeros(len(vect)) #Rényi relative entropy
TRE = np.zeros(len(vect)) #Tsallis relative entropy
# R = np.zeros(len(vect)) #Rate function


alfa = 0.8


for tt in range(len(vect)):
    U = expm(-1j*Hf*vect[tt])  #time evolution operator
    Ut = np.transpose(np.conjugate(U))
    
    rhot = np.dot(U, np.dot(rho0, Ut))  #quenched density matrix at time tt
    
    # RE[tt] = np.trace(np.dot(rhot, logm(rhot) - logm(rhof_eq)))
    # RE[tt] = np.trace(np.dot(rhot, beta*(Hf-H0))) + np.log(Zf/Z0)
    # RRE[tt] = 1/(alfa-1)*np.log(np.trace(np.dot(fractional_matrix_power(rhot, alfa), fractional_matrix_power(rho0, 1-alfa))))
    TRE[tt] = 1/(1-alfa)*(1-np.trace(np.dot(fractional_matrix_power(rhot, alfa), fractional_matrix_power(rho0, 1-alfa))))
    


# #find_peaks command helps to find the peaks of the derivative of relative function
# # tm, _ = find_peaks(np.diff(TRE)/(vect[1]-vect[0]), height=-150)


# # plt.plot(gx*(vect[1:]+vect[0:-1])/2, np.diff(TRE)/(vect[1]-vect[0]), label=r"$h={}, \beta={}, \alpha={}$".format(h, beta, alfa))
# # plt.plot(gx*(vect[1:]+vect[0:-1])[tm]/2, np.diff(TRE)[tm]/(vect[1]-vect[0]), 'x')

plt.plot(gx*vect, TRE, label=r"$\gamma={}, h={}, \beta={}, j={}$".format(gx, h, beta, j))
plt.xlabel("$\gamma_x t$", fontsize=17)
# plt.xlim([0,30])
plt.ylabel(r"$R(\rho(t)||\rho_0)$", fontsize=17)
plt.legend(fontsize=12)
plt.tight_layout()


# plt.savefig("RenyiEntropy_ThermalState_j100_iso.pdf")


# plt.axvline(8.02675585284281, c='r', linewidth=1) j=100
# plt.axvline(12.94314381270903, c='r', linewidth=1)
# plt.axvline(17.959866220735787, c='r', linewidth=1)
# plt.axvline(22.876254180602007, c='r', linewidth=1)
# plt.axvline(8.02675585284281, c='r', linewidth=1) #j=300
# plt.axvline(12.94314381270903, c='r', linewidth=1)
# plt.axvline(17.959866220735787, c='r', linewidth=1)
# plt.axvline(22.876254180602007, c='r', linewidth=1)