
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 100
dim = np.int(2*j+1)
vect = np.linspace(0,50,350)



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

def expval(v1,A,v2):
    v1_dagger = np.transpose(np.conjugate(v1))
    return np.dot(v1_dagger, np.dot(A, v2))

#Initial hamiltonian

omega = 0
chi=1
w=0
H0 = -(chi/(2*j))*np.dot(Jz(j),Jz(j)) - omega*Jx(j) - w*Jz(j)


#Initial density matrix of the system
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index]
V0 = V0[:,index]


# Initial state
psi0 = (1/np.sqrt(2))*(V0[:,0] + V0[:,1])
# psi0 = V0[:,0]


#Final hamiltonian
omega = 0.3
chi=1
w=0
H = -(chi/(2*j))*np.dot(Jz(j),Jz(j)) - omega*Jx(j) - w*Jz(j)


#####################################################################
#
# Loschmidt Echo
#
#####################################################################

L = np.zeros(len(vect))  #Loschmidt Echo
# Precision = np.zeros(len(vect))

for tt in range(len(vect)):
    U = expm(-1j*H*vect[tt])
    
    L[tt] = np.abs(expval(psi0, U, psi0))**2
    
    
    # ################ precision ##############
    # psit = np.dot(U, V0[:,0])
    # F = Jz(j)
    
    # Mean0 = np.dot(dag(V0[:,0]), np.dot(F, V0[:,0]))
    # Var0 = np.sqrt(np.dot(dag(V0[:,0]), np.dot(np.dot(F,F), V0[:,0])) - Mean0**2)
    
    # Mean_t = np.dot(dag(psit), np.dot(F, psit))
    # Var_t = np.sqrt(np.dot(dag(psit), np.dot(np.dot(F,F), psit)) - Mean_t**2)
    
    # Precision[tt] = ((Var0 + Var_t)/(Mean0 - Mean_t))**2
    # # R[tt] = -1/(2*j)*np.log(L)
    


# Plotting
plt.plot(vect, L, label=r"$N={}$".format(2*j))
plt.xlabel(r"$t$", fontsize=17)
plt.ylabel(r"$\mathcal{L}(t)$", fontsize=17)
plt.legend(fontsize=12)
plt.tight_layout()





