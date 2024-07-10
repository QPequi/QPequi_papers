


#######################################
#
# This code calculates the Krylov complexity for
# a quench dynamics in the LMG model
#
#######################################

import numpy as np
from numpy import linalg as LA
from scipy.linalg import expm
from math import isnan
# from functools import reduce
# from scipy.signal import find_peaks
import matplotlib.pyplot as plt


import warnings
warnings.filterwarnings('ignore')

plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 1000
dim = int(2*j+1)
g = 0.1



#Defining the Jx, Jy and Jz in z basis and some other useful functions
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

            p1 = int(m1+j)
            p2 = int(m2+j)

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

            p2 = int(m2+j)
            p1 = int(m1+j)

            Jy[p1,p2] = 1j*(aut1-aut2)

            aut1 = aut2 = 0
    return Jy

def Jz(j):
    Jz = np.zeros((dim, dim))

    for m in np.arange(-j, j+1, 1):
        p = int(m+j)
        Jz[p,p] = -m
    return Jz

# Commutator
def Comm(A,B):
    return np.dot(A,B)-np.dot(B,A)

#Bra-ket 
def mean(v1,A,v2):
    v1_dagger = np.transpose(np.conjugate(v1))
    return np.dot(v1_dagger, np.dot(A, v2))

    
###############################
#
# Quenching the systems
#
###############################
#Initial hamiltonian
h0 = 0
H0 = -1/j*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h0*Jz(j)


#Initial density matrix of the system
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index]
V0 = V0[:,index]

#Final hamiltonian
h = 0.3
H = -1/j*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h*Jz(j)


################################
#
# Generating the Krylov basis
#
################################

N = 50

a = np.zeros(N) #Lanczos coefficients
b = np.zeros(N)

K = np.zeros((dim,N))  # K stores the Krylov states in the columns

K[:,0] = V0[:,0]  #Lanczos algorithm initiates with the initial ground state
a[0] = mean(V0[:,0], H, V0[:,0])

K[:,1] = np.dot(H, V0[:,0]) - a[0]*V0[:,0]
b[1] = LA.norm(K[:,1])  #normalizing K_1
K[:,1] = K[:,1]/b[1]

for j in range(2, N):  #Producing the next Krylov states
    a[j-1] = np.dot(K[:,j-1], np.dot(H, K[:,j-1]))
    K[:,j] = np.dot(H, K[:,j-1]) - a[j-1]*K[:,j-1] - b[j-1]*K[:,j-2] #Lanczos algorithm
    
    b[j] = LA.norm(K[:,j]) #Normalization
    K[:,j] = K[:,j]/b[j]


##################################
#
# Calculating the complexity
#
##################################

vect = np.linspace(0, 30, 200)

C = np.zeros(len(vect)) #Stores the complexity

for tt in range(len(vect)):
    U = expm(-1j*H*vect[tt])  #Time evolution operator

    psit = np.dot(U, V0[:,0])

    for n in range(N):  #Complexity

        C[tt] = C[tt] + n*np.abs(np.dot(np.conjugate(psit), K[:,n]))**2

#Plotting
plt.plot(vect, C, label=r"$j={}$".format(j))
plt.xlabel(r"$t$", fontsize=17)
plt.ylabel(r"$C(t)$", fontsize=17)
plt.tight_layout()





