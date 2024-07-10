

import numpy as np
from numpy import linalg as LA
from functools import reduce
from scipy.linalg import expm
from scipy.linalg import fractional_matrix_power
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


N = 10
J = 1
beta = 1

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
h = 0

H0 = 0
for j in range(N-1):  #the first site is j=0 and the last j=N-1
    H0 = H0 -J*np.dot(sx(j), sx(j+1)) - h*sz(j) 
    
H0 = H0 -J*np.dot(sx(N-1), sx(0)) - h*sz(N-1)  #periodic boundary conditions


#Diagonalizing the initial hamiltonian
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index] 
V0 = V0[:,index]



#Final Ising hamiiltonian
h = 1.2

Hf = 0
for j in range(N-1):  #the first site is j=0 and the last j=N-1
    Hf = Hf -J*np.dot(sx(j), sx(j+1)) - h*sz(j) 
    
Hf = Hf -J*np.dot(sx(N-1), sx(0)) - h*sz(N-1)  #periodic boundary conditions


#########################################
#
# Loschmidt echo and rate function
#
#########################################

vect = np.linspace(0, 30, 200)

LE = np.zeros(len(vect))

for tt in range(len(vect)):
    
    U = expm(-1j*vect[tt]*Hf)
    
    LE[tt] = np.abs(np.dot(V0[:,0], np.dot(U, V0[:,0])))**2
    

plt.plot(vect, LE)

