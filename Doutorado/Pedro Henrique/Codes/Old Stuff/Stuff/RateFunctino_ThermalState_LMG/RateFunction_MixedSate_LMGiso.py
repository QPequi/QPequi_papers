


import numpy as np
from scipy.linalg import expm
from scipy import linalg as LA
import matplotlib.pyplot as plt
from scipy.signal import find_peaks #it helps to identify t_c (peaks)


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 100
dim = np.int(2*j+1)
gx = 1
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


#Commutator [A,B]
def Comm(A,B):
    return np.dot(A,B)-np.dot(B,A)


# Energies of LMG hamiltonian
def E(J, mz, h):
    return -2*h*mz - (gx/J)*(J*(J+1) - mz**2)


#Initial LMG hamiltonian
h = 0
H0 = -2*h*Jz(j) - (gx/j)*(np.dot(Jx(j),Jx(j)) + np.dot(Jy(j),Jy(j)))


# Initial hamiltonian's eigenvalues and eigenvectors
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index]
V0 = V0[:,index]


beta = 2


Z0 = 0
for m in range(-j, j+1, 1):
    Z0 = Z0 + np.exp(-beta*E(j,m,h))


# standard basis vectors
Vz = np.zeros((dim,dim))
for i in range(dim):
    Vz[dim-1-i, i] = 1
    

# # # # Initial density matrix: ground state/thermal state
# rho0 = np.dot(V0[:,0].reshape(dim,1), V0[:,0].reshape(1,dim))
rho0 = np.zeros((dim, dim))
for n in range(-j, j+1, 1):
    rho0 = rho0 + (1/Z0)*np.exp(-beta*E(j, n, h))*np.dot(Vz[:,n+j].reshape(dim,1), Vz[:,n+j].reshape(1,dim))


h = 0.8

#Final LMG hamiltonian
Hf = -2*h*Jz(j) - (gx/j)*(np.dot(Jx(j),Jx(j)) + np.dot(Jy(j),Jy(j)))

# Zf = 0
# for m in range(-j, j+1, 1):
#     Zf = Zf + np.exp(-beta*E(j,m,gx))


# rhof_eq = np.zeros((dim, dim))
# for n in range(-j, j+1, 1):
#     rhof_eq = rhof_eq + (1/Zf)*np.exp(-beta*E(j, n, gx))*np.dot(Vz[:,n+j].reshape(dim,1), Vz[:,n+j].reshape(1,dim))




######################################################
#
# Rate function for thermal initial state
#
######################################################


LE_n = np.zeros(len(vect)) #Loschmidt echo
LE_a = np.zeros(len(vect)) #Loschmidt echo
# R = np.zeros(len(vect)) #rate function


for tt in range(len(vect)):
    U = expm(-1j*vect[tt]*Hf)
    
    LE_n[tt] = np.abs(np.trace(np.dot(rho0, U)))**2
    # R[tt] = -1/(2*j)*np.log(LE)
    
    # # Initial pure state:
    # O = 0
    # for m in range(-j, j+1, 1):
    #     O = O + np.exp(-1j*E(j, m, h)*vect[tt])*np.abs(np.dot(Vz[:,m+j], Vz[:,0]))**2
    
    # LE_a[tt] = np.abs(O)**2
    
    # Initial mixed state:
    # O = 0
    # for m in range(-j, j+1, 1):
    #     # for n in range(-j, j+1, 1):
    #     O = O + (1/Z0)*np.exp(-beta*E(j, m, 0))*np.exp(-1j*E(j, m, h)*vect[tt])

    # LE_a[tt] = np.abs(O)**2

#find_peaks command helps to identify where the peaks occur.
# tc, _ = find_peaks(R, height=0)


plt.plot(vect, LE_n, label=r"$j={}$".format(j))
# plt.plot(vect, LE_a, label=r"$j={}$".format(j))
# plt.plot(vect, -1/(2*j)*np.log(LE_a), label=r"$j={}$".format(j))
plt.xlabel("$\gamma_x t$", fontsize=17)
plt.ylabel(r"$\mathcal{L}(t)$", fontsize=17)
plt.legend(fontsize=10, loc=1)
plt.tight_layout()
    

plt.savefig("Analyt_LoschmidtEcho_PureInitState_h0.8.pdf")










