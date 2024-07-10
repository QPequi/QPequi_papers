



import numpy as np
from scipy.linalg import expm
from scipy import linalg as LA
from scipy.linalg import fractional_matrix_power
import matplotlib.pyplot as plt
# from scipy.signal import find_peaks #it helps to identify t_c (peaks)


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 100
dim = np.int(2*j+1)
g = 1
vect = np.linspace(0, 30, 300)



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

# Energies of LMG hamiltonian
def E(J, jz, h):
    return -2*h*jz - (g/J)*(J*(J+1) - jz**2) + g

#Commutator 
def Comm(A,B):
    return np.dot(A,B)-np.dot(B,A)


#Schatten 2-norm
def norm2(A):
    return np.sqrt(np.trace(np.dot(np.transpose(np.conjugate(A)), A)))


#Initial LMG hamiltonian
h = 0
H0 = -2*h*Jz(j) - (g/j)*(np.dot(Jx(j),Jx(j)) + np.dot(Jy(j),Jy(j))) + g


#Initial hamiltonian's engeinvalues and eigenvectors
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index]
V0 = V0[:,index]


beta = 2

Z0 = 0

for m in range(-j, j+1, 1):
    Z0 = Z0 + np.exp(-beta*E(j,m,h))


# Diagonalizing Jz to construct Dicke states
Ez, Vz = LA.eig(Jz(j))
index = np.argsort(Ez) #Get the index list that would sort E and V
Ez = Ez[index]
Vz = Vz[:,index]


# # Initial density matrix: ground state/thermal state
# rho0 = np.dot(V0[:,0].reshape(dim,1), V0[:,0].reshape(1,dim))
rho0 = np.zeros((dim, dim))
for n in range(-j, j+1, 1):
    rho0 = rho0 + (1/Z0)*np.exp(-beta*E(j, n, h))*np.dot(Vz[:,n+j].reshape(dim,1), Vz[:,n+j].reshape(1,dim))


h = 2

#Final LMG hamiltonian
Hf = -2*h*Jz(j) - (g/j)*(np.dot(Jx(j),Jx(j)) + np.dot(Jy(j),Jy(j))) + g

# Zf = 0
# for m in range(-j, j+1, 1):
#     Zf = Zf + np.exp(-beta*E(j,m,gx))



# rhof_eq = np.zeros((dim, dim))
# for n in range(-j, j+1, 1):
#     rhof_eq = rhof_eq + (1/Zf)*np.exp(-beta*E(j, n, gx))*np.dot(Vz[:,n+j].reshape(dim,1), Vz[:,n+j].reshape(1,dim))


#Final density matrix of the system
Ef, Vf = LA.eig(Hf)
index = np.argsort(Ef) #Get the index list that would sort E and V
Ef = Ef[index]
Vf = Vf[:,index]

######################################################
#
# Rate function for thermal initial state
#
######################################################


# RRE = np.zeros(len(vect)) #Rényi relative entropy
# TRE = np.zeros(len(vect)) #Tsallis relative entropy

# tau = np.zeros(len(vect)) #QSL time

# alfa = 0.4

# G = np.abs(1+(1-alfa)*np.log(np.min(E0)))**(-1)*norm2(fractional_matrix_power(rho0, 1-alfa))*norm2(Comm(Hf, fractional_matrix_power(rho0, alfa)))
# # G = norm2(fractional_matrix_power(rho0, 1-alfa))*norm2(Comm(Hf, fractional_matrix_power(rho0, alfa)))
    

# for tt in range(len(vect)):
#     U = expm(-1j*Hf*vect[tt])  #time evolution operator
#     Ut = np.transpose(np.conjugate(U))
    
#     rhot = np.dot(U, np.dot(rho0, Ut))  #quenched density matrix at time tt
    
#     RRE = 1/(alfa-1)*np.log(np.trace(np.dot(fractional_matrix_power(rhot, alfa), fractional_matrix_power(rho0, 1-alfa))))
#     # TRE = 1/(1-alfa)*(1-np.trace(np.dot(fractional_matrix_power(rhot, alfa), fractional_matrix_power(rho0, 1-alfa))))
    

#     tau[tt] = np.abs(1-alfa)*RRE/G


# #find_peaks command helps to find the peaks of the derivative of relative function
# # tc, _ = find_peaks(np.diff(TRE)/(vect[1]-vect[0]), height=-150)


# # plt.plot(gx*(vect[1:]+vect[0:-1])/2, np.diff(TRE)/(vect[1]-vect[0]), label=r"$h={}, \beta={}, \alpha={}$".format(h, beta, alfa))
# # plt.plot(gx*(vect[1:]+vect[0:-1])[tm]/2, np.diff(TRE)[tm]/(vect[1]-vect[0]), 'x')

# plt.plot(vect, tau, label=r"$j={}$".format(j))
# plt.xlabel("$\gamma_x t$", fontsize=17)
# # plt.xlim([0,30])
# plt.ylabel(r"$\tau^{R}_{\alpha}(\rho(t)||\rho_0)$", fontsize=17)
# plt.legend(fontsize=12)
# plt.tight_layout()

#np.savetxt("RateFunction_j50_beta1_h2.txt", R)









