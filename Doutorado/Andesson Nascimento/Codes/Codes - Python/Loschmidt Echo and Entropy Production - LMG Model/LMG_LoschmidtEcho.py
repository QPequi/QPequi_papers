import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
import math as mt


# plt.rcParams["font.size"] = 15
# plt.rcParams["text.usetex"] = True

j = 1000
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
h = 0.8
H = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h*Jz(j)


#####################################################################
#
# Loschmidt Echo
#
#####################################################################

L = np.zeros(len(vect))  #Loschmidt Echo
LL = np.zeros(len(vect)) #Bures Angle
EP = np.zeros(len(vect)) #Minimal Entropy Production
EP1 = np.zeros(len(vect))
EP2 = np.zeros(len(vect))
EP3 = np.zeros(len(vect))
EP4 = np.zeros(len(vect))
RF = np.zeros(len(vect)) #Rate Function
DEP = np.zeros(len(vect)) #Variation Entropy Production

for tt in range(len(vect)):
    U = expm(-1j*H*vect[tt])
    
    L[tt] = np.abs(np.dot(dag(V0[:,0]), np.dot(U, V0[:,0])))**2
    LL[tt] = mt.acos(L[tt])
    EP[tt] = 8/(mt.pi**2)*LL[tt]**2 + 64/(9*mt.pi**4)*LL[tt]**4 + \
    (32*(2**6))/(135*mt.pi**6)*LL[tt]**6 + (992*(2**8))/(5103*mt.pi**8)*LL[tt]**8 + \
    (6656*(2**10))/(32805*mt.pi**10)*LL[tt]**10
    EP4[tt] = 8/(mt.pi**2)*LL[tt]**2 + 64/(9*mt.pi**4)*LL[tt]**4 + \
    (32*(2**6))/(135*mt.pi**6)*LL[tt]**6 + (992*(2**8))/(5103*mt.pi**8)*LL[tt]**8
    EP3[tt] = 8/(mt.pi**2)*LL[tt]**2 + 64/(9*mt.pi**4)*LL[tt]**4 + \
    (32*(2**6))/(135*mt.pi**6)*LL[tt]**6 
    EP2[tt] = 8/(mt.pi**2)*LL[tt]**2 + 64/(9*mt.pi**4)*LL[tt]**4
    EP1[tt] = 8/(mt.pi**2)*LL[tt]**2
    RF[tt] = -(1/(2*j))*mt.log(L[tt])
 
dEP = np.gradient(EP)
dt = np.gradient(vect)   
# dEP = np.diff(EP)
# dt = np.diff(vect)

DEP = dEP/dt



# Plotting
# plt.plot(vect, L, label="Loschmidt Echo")
# #plt.plot(vect, LL, label=r"$N={}$".format(2*j))
# ##plt.plot(vect, EP, label="Entropy Production")
# #plt.plot(vect, EP, label=r"$N={}$".format(2*j))
# plt.xlabel(r"$t$", fontsize=17)
# #plt.ylabel(r"$\mathcal{L}(t)$", fontsize=17)
# plt.legend(fontsize=12)
# plt.tight_layout()

fig, ax = plt.subplots(2,2,figsize =(15,10))
ax[0][0].plot(vect,L,label='Loschmidt Echo') # 1º grafico
# ax[0][0].ylabel("L",fontsize=17)
# ax[0][0].xlabel(r"$t$", fontsize=17)
ax[0][0].legend()
ax[0][0].grid()
ax[0][1].plot(vect,RF,label='Rate Function') # 2º grafico
# ax[0][1].ylabel("RF",fontsize=17)
# ax[0][1].xlabel(r"$t$", fontsize=17)
ax[0][1].legend()
ax[0][1].grid()
ax[1][0].plot(vect,EP,label='Entropy Production') # 3º grafico
ax[1][0].plot(vect,EP1,label='Entropy Production1')
ax[1][0].plot(vect,EP2,label='Entropy Production2')
ax[1][0].plot(vect,EP3,label='Entropy Production3')
ax[1][0].plot(vect,EP4,label='Entropy Production4')
# ax[1][0].ylabel("EP",fontsize=17)
# ax[1][0].xlabel(r"$t$", fontsize=17)
ax[1][0].legend()
ax[1][0].grid()
ax[1][1].plot(vect,DEP,label='Variation of EP') # 4º grafico
# ax[1][1].ylabel("dS/dt",fontsize=17)
# ax[1][1].xlabel(r"$t$", fontsize=17)
ax[1][1].legend()
ax[1][1].grid()



