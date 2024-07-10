


import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 100
dim = np.int(2*j+1)

vect = np.linspace(0, 30, 150)



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


#Initial hamiltonian
h = 1
H0 = -(1-h)*Jz(j) - h/2/j*np.dot(Jx(j),Jx(j))


#Initial hamiltonian's eigenvalues and eigenvectors
E, V = LA.eig(H0)
index = np.argsort(E) #Get the index list that would sort E and V
E = E[index] 
V = V[:,index]


#initial state
psi0 = V[:,0]  #pure ground-state of H0


#Final hamiltonian
h = 0.5
Hf = -(1-h)*Jz(j) - h/2/j*np.dot(Jx(j),Jx(j))


#####################################################################
#
# Loschmidt Echo
#
#####################################################################

L = np.zeros(len(vect))  #Loschmidt Echo

for tt in range(len(vect)):
    U = expm(-1j*Hf*vect[tt])
    
    L[tt] = np.abs(np.dot(psi0, np.dot(U, psi0)))**2
    
    # R[tt] = -1/(2*j)*np.log(L)


plt.plot(vect, L, label='$h={}, j={}$'.format(h,j))
# plt.plot(vect[tc], L[tm], 'x')
plt.legend(fontsize=12)
plt.xlabel('$\gamma_xt$', fontsize=18)
plt.ylabel('$\mathcal{L}(t)$', fontsize=18)
plt.tight_layout()

# plt.savefig('LE_h0.5_j100.pdf')



#np.savetxt("LE_j200.txt", psi0)





# plt.plot(gx*vect, tau, label=r"$\tau^{R}_{\alpha}(\rho(t)||\rho_0)$")
# plt.plot(gx*vect, 20*L, label=r"$\mathcal{L}(t)$")
# plt.xlabel(r"$\gamma_x t$", fontsize=17)
# plt.title(r"$h={},\alpha={},j={}$".format(h,alfa,j))
# plt.legend(fontsize=12)
# plt.tight_layout()
# plt.savefig("QSL_LE_Renyi_purestate_j300.pdf")

