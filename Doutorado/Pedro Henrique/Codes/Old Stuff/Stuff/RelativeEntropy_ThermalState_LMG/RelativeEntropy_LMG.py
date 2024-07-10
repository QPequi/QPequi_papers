


import numpy as np
from numpy import linalg as LA
from scipy.linalg import fractional_matrix_power
from scipy.linalg import expm
from scipy.linalg import logm
from scipy.signal import find_peaks
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 100
dim = np.int(2*j+1)
gx = 1
gy = 0
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

#Commutator 
def Comm(A,B):
    return np.dot(A,B)-np.dot(B,A)
    

#Initial LMG hamiltonian
h=0
H0 = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))


#Initial density matrix of the system
E, V = LA.eig(H0)
index = np.argsort(E) #Get the index list that would sort E and V
E = E[index] 
V = V[:,index]

beta = 1

#Initial density matrix: ground state/thermal state
# rho0 = np.dot(V[:,0].reshape(dim, 1), V[:,0].reshape(1, dim))
rho0 = expm(-beta*H0)/np.trace(expm(-beta*H0))


h = 1.2

#Final LMG hamiltonian
Hf = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))


#Final equilibrium state 
rhof_eq = expm(-beta*Hf)/np.trace(expm(-beta*Hf))


####################################################
#
# Relative entropy 
#
####################################################


RE = np.zeros(len(vect)) #usual relative entropy
# RRE = np.zeros(len(vect)) #Rényi relative entropy
# TRE = np.zeros(len(vect)) #Tsallis relative entropy
# R = np.zeros(len(vect)) #Rate function


alfa = 0.3

for tt in range(len(vect)):
    U = expm(-1j*Hf*vect[tt])  #time evolution operator
    Ut = np.transpose(np.conjugate(U))
    
    rhot = np.dot(U, np.dot(rho0, Ut))  #quenched density matrix at time tt
    
    RE[tt] = np.trace(np.dot(rhot, logm(rhot) - logm(rhof_eq)))
    # RE[tt] = np.trace(np.dot(rhot, beta*(Hf-H0))) + np.log(np.trace(expm(-beta*Hf))/np.trace(expm(-beta*H0)))
    # RRE[tt] = 1/(alfa-1)*np.log(np.trace(np.dot(fractional_matrix_power(rhot, alfa), fractional_matrix_power(rhof_eq, 1-alfa))))
    # TRE[tt] = 1/(1-alfa)*(1-np.trace(np.dot(fractional_matrix_power(rhot, alfa), fractional_matrix_power(rhof_eq, 1-alfa))))
    



#find_peaks command helps to find the peaks of the derivative of relative function
# tm, _ = find_peaks(np.diff(TRE)/(vect[1]-vect[0]), height=-150)


# plt.plot(gx*(vect[1:]+vect[0:-1])/2, np.diff(TRE)/(vect[1]-vect[0]), label=r"$h={}, \beta={}, \alpha={}$".format(h, beta, alfa))
# plt.plot(gx*(vect[1:]+vect[0:-1])[tm]/2, np.diff(TRE)[tm]/(vect[1]-vect[0]), 'x')

plt.plot(gx*vect, RE, label=r"$h={}, \beta={}, \alpha={}, j={}$".format(h, alfa, beta, j))
# plt.plot(gx*vect, TRE[tm], 'x')
plt.xlabel("$\gamma_x t$", fontsize=17)
# plt.xlim([0,30])
plt.ylabel(r"$H(\rho(t)||\rho_f^{eq})$", fontsize=17)
plt.legend(loc=4, fontsize=12)
plt.tight_layout()


# plt.savefig("TsallisEntropy_ThermalState_j100_alpha0.3.pdf")


# plt.axvline(8.02675585284281, c='r', linewidth=1) j=100
# plt.axvline(12.94314381270903, c='r', linewidth=1)
# plt.axvline(17.959866220735787, c='r', linewidth=1)
# plt.axvline(22.876254180602007, c='r', linewidth=1)
# plt.axvline(8.02675585284281, c='r', linewidth=1) #j=300
# plt.axvline(12.94314381270903, c='r', linewidth=1)
# plt.axvline(17.959866220735787, c='r', linewidth=1)
# plt.axvline(22.876254180602007, c='r', linewidth=1)


# plt.subplot(211)
# plt.plot(gx*vect, RRE)
# plt.ylabel(r"$R(\rho(t)||\rho_f)$", fontsize=17)
# plt.xlabel(r"$\gamma_x t$", fontsize=17)
# plt.tight_layout()

# plt.subplot(212)
# plt.plot(gx*vect, tau)
# plt.ylabel(r"$\tau^{R}_{\alpha}(\rho(t)||\rho_f)$", fontsize=17)
# plt.xlabel(r"$\gamma_x t$", fontsize=17)
# plt.tight_layout()

# plt.savefig("QSLtime_RelatEntro_Renyi_j100.pdf")
# label=r"$\alpha={}, \beta={}, h={}, j={}$".format(alfa,beta,h,j)