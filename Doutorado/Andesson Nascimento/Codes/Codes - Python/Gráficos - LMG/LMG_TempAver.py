import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
import math as mt


# plt.rcParams["font.size"] = 15
# plt.rcParams["text.usetex"] = True

ll = 1
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

#Final hamiltonian
h = 0.2
H = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h*Jz(j)

#####################################################################
#
# Loschmidt Echo
#
#####################################################################

L = np.zeros(len(vect))  # Loschmidt Echo
LL = np.zeros(len(vect)) # Bures Angle
EP = np.zeros(len(vect)) # Minimal Entropy Production - 1 term
RFN = np.zeros(len(vect)) # Rate Function * ll ("Normalized")
RF = np.zeros(len(vect)) # Rate Function
DEP = np.zeros(len(vect)) # Variation Entropy Production
DEPN = np.zeros(len(vect)) # Variation Entropy Production * 0,1 ("Normalized")
TAEP = np.zeros(len(vect)) #  Time Average - Entropy production

for tt in range(len(vect)):
    U = expm(-1j*H*vect[tt])
    
    L[tt] = np.abs(np.dot(dag(V0[:,0]), np.dot(U, V0[:,0])))**2 
    LL[tt] = mt.acos(L[tt])
    (6656*(2**10))/(32805*mt.pi**10)*LL[tt]**10
    EP[tt] = 8/(mt.pi**2)*LL[tt]**2
    RFN[tt] = -(1/(2*j))*mt.log(L[tt])*ll
    RF[tt] = -(1/(2*j))*mt.log(L[tt])
    TAEP[tt] = (1/np.max(vect))*np.trapz(EP, vect)

# for TT in range(len(vect)):
    # TAEP[TT] = (1/np.max(vect))*np.trapz(EP[tt], vect)

# Variation of Entropy Production

dEP = np.gradient(EP)
dt = np.gradient(vect)   

DEP = dEP/dt
DEPN = (dEP/dt)*0.1

# Plotting

fig, ax = plt.subplots(2,2,figsize =(15,10))
ax[0][0].plot(vect,L,label='Loschmidt Echo') # 1º grafico
ax[0][0].plot(vect,RFN, label='Rate Function') # 1º grafico
ax[0][0].plot(vect,EP,label='Entropy Production') # 2º grafico
ax[0][0].legend()
ax[0][0].grid()

ax[0][1].plot(vect,L,label='Loschmidt Echo') # 2º grafico
ax[0][1].plot(vect,LL,label='Bures Angle ') # 3º grafico
ax[0][1].legend()
ax[0][1].grid()

ax[1][0].plot(vect,RF,label='Rate Function') # 3º grafico
ax[1][0].plot(vect,DEPN,label='Variation of EP(x0,1)') # 3º grafico
ax[1][0].legend()
ax[1][0].grid()

ax[1][1].plot(vect,TAEP,label='Temporal Average of Entropy Production') # 4º grafico
ax[1][1].legend()
ax[1][1].grid()



