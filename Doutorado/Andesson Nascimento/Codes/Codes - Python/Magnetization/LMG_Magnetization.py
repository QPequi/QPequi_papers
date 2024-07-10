import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
import math as mt


# plt.rcParams["font.size"] = 15
# plt.rcParams["text.usetex"] = True

ll = 2
j = ll*100
dim = np.int(2*j+1)
g = 0.1
vect = np.linspace(0,30,200)

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

h0 = 0
H0 = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h0*Jz(j)

#Initial density matrix of the system
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index]
V0 = V0[:,index]

# #Final hamiltonian
# h = 0.6
# H = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h*Jz(j)

#####################################################################
#
# Loschmidt Echo
#
#####################################################################

L = np.zeros(len(vect))  # Loschmidt Echo
RF = np.zeros(len(vect)) # Rate Function
TAEP = np.zeros(len(vect)) #  Time Average - Entropy production
M = np.zeros(len(vect)) #magnetization vector

# vech = np.linspace(0.4,0.6,2)
# for h in range(len(vech)):
#     H = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h*Jz(j)
#     for tt in range(len(vect)):
#         U = expm(-1j*H*vect[tt])
#         psit = np.dot(U, V0[:,0])
#         L[tt] = np.abs(np.dot(dag(V0[:,0]), np.dot(U, V0[:,0])))**2 
#         RF[tt] = -(1/(2*j))*mt.log(L[tt])
#         M[tt] = np.dot(np.conjugate(psit), np.dot(Jx(j), psit)) #magnetization at time vect[tt]
#     # TAEP[tt] = (1/np.max(vect))*np.trapz(EP, vect)
#     plt.plot(vect, M, label=h)

h=0.8
H = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h*Jz(j)
for tt in range(len(vect)):
    U = expm(-1j*H*vect[tt])
    psit = np.dot(U, V0[:,0])
    L[tt] = np.abs(np.dot(dag(V0[:,0]), np.dot(U, V0[:,0])))**2 
    RF[tt] = -(1/(2*j))*mt.log(L[tt])
    M[tt] = np.dot(np.conjugate(psit), np.dot(Jx(j), psit)) #magnetization at time vect[tt]
    


plt.plot(vect, M)









# Plotting


# plt.plot(vect, function, label=r"name")
# plt.legend(fontsize=15)
# plt.xlabel('$\gamma_xt$', fontsize=18)
# plt.ylabel('$\mathcal{L}(t)$', fontsize=18)
# plt.tight_layout()

# plt.xlim([-0.2, 20])

# plt.savefig('LE_h0.5_j100.pdf')

# fig, ax = plt.subplots()
# left, bottom, width, height = [.30, 0.6, 0.2, 0.25]
# ax_new = fig.add_axes([left, bottom, width, height])
# ax.plot(x, y, color='red')
# ax_new.plot(x, y, color='green')
# plt.show()