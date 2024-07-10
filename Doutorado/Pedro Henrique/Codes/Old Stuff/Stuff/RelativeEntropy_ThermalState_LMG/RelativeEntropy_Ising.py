

import numpy as np
from functools import reduce
from scipy.linalg import expm
from numpy import linalg as LA
from scipy.linalg import fractional_matrix_power
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


N = 2
J = 1
beta = 0.1

# sigmaZ_j: PauliMatrix(z) of site j
def sz(j):
    S = []
    for i in range(N): 
        S.append(np.eye(2))
    
    S[j] = np.array([[1,0],[0,-1]])
    S = reduce(np.kron, S)
    
    return S

# sigmaX_j: PauliMatrix(x) of site j
def sx(j): 
    S = []
    for i in range(N):
        S.append(np.eye(2))
    
    S[j] = np.array([[0,1],[1,0]])
    S = reduce(np.kron, S)
    
    return S

#Initial Ising hamiltonian
h = 0.2

H0 = 0
for j in range(N-1):  #the first site is j=0 and the last j=N-1
    H0 = H0 -J*np.dot(sx(j), sx(j+1)) - h*sz(j) 
    
H0 = H0 -J*np.dot(sx(N-1), -sx(0)) - h*sz(N-1)  #(anti)periodic boundary conditions



#Initial density matrix of the system
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index] 
V0 = V0[:,index]


#Initial thermal state
rho0 = expm(-beta*H0)/np.trace(expm(-beta*H0))
# rho0 = np.dot(V0[:,0].reshape(2**N, 1), V0[:,0].reshape(1, 2**N))


#Final Ising hamiiltonian
h = 0.8

Hf = 0
for j in range(N-1):  #the first site is j=0 and the last j=N-1
    Hf = Hf -J*np.dot(sx(j), sx(j+1)) - h*sz(j) 
    
Hf = Hf -J*np.dot(sx(N-1), sx(0)) - h*sz(N-1)  #periodic boundary conditions


    
##############################################
#
# Renyi/Tsallis relative entropy
#
##############################################

vect = np.linspace(0, 30, 200)

RRE = np.zeros(len(vect)) #Rényi relative entropy

alpha = 0.8

for tt in range(len(vect)):
    
    U = expm(-1j*vect[tt]*Hf)  #time evolution operator
    Ud = np.transpose(np.conjugate(U))  #U^\dagger
    
    rhot = np.dot(U, np.dot(rho0, Ud))  #rho(t)
    
    RRE[tt] = 1/(alpha-1)*np.log(np.trace(np.dot(fractional_matrix_power(rhot, alpha), fractional_matrix_power(rho0, 1-alpha))))
    

plt.plot(vect, RRE, label=r"$N={}, \beta={}, \alpha={}$".format(N,beta,alpha))
# plt.plot(vect, LE, label=r"$LE$")
plt.xlabel("$t$", fontsize=17)
plt.ylabel(r"$R(\rho(t)||\rho_0)$", fontsize=17)
plt.title("Ising\ Model")
plt.legend(fontsize=10)
plt.tight_layout()

# plt.savefig("RenyiRelativeEntropy_IsingModel.pdf")















