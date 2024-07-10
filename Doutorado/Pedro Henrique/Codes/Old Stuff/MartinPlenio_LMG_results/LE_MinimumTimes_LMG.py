

#########################################
#
# Testing the expression of the LE minimal times
# when lambda_i, lambda_f < 1 (equation S16 from arXiv:1904.09937)
#
#########################################



import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 300
dim = np.int(2*j+1)
vect = np.linspace(0, 50, 150)


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


########## Fixed parameters ########

J = 1
g = 1

####################################

#Initial hamiltonian
h0 = 1.2
H0 = -2*J/(2*j)*(np.dot(Jx(j), Jx(j)) + g*np.dot(Jy(j), Jy(j))) - 2*h0*Jz(j) + (1+g)*J/2*np.identity(dim)


#Initial hamiltonian's eigenvalues and eigenvectors
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index] 
V0 = V0[:,index]


#initial state
psi0 = V0[:,0]  #pure ground-state of H0

#Final hamiltonian
h = 1.3
Hf = -2*J/(2*j)*(np.dot(Jx(j), Jx(j)) + g*np.dot(Jy(j), Jy(j))) - 2*h*Jz(j) + (1+g)*J/2*np.identity(dim)



###################################################
#
# Loschmidt echo
#
###################################################

LE = np.zeros(len(vect))  #Loschmidt Echo


for tt in range(len(vect)):
    U = expm(-1j*Hf*vect[tt])
    
    psit = np.dot(U, psi0)
    
    LE[tt] = np.abs(np.dot(np.conjugate(psi0), psit))**2



plt.plot(vect, LE, label=r"$N={}, h_i={}, h_f={}$".format((2*j), h0, h))
plt.legend(fontsize=15)
plt.title(r"$J={}, \gamma={}$".format(J, g))
plt.xlabel(r'$t$', fontsize=18)
# plt.ylabel(r'$\overline{S_d}$', fontsize=18)
plt.tight_layout()
# plt.savefig("TimeAverage_DiagEntropy_N1000_ParametLMG.pdf")

# plt.savefig('LE_h0.5_j100.pdf')



#np.savetxt("LE_j200.txt", psi0)





# plt.plot(gx*vect, tau, label=r"$\tau^{R}_{\alpha}(\rho(t)||\rho_0)$")
# plt.plot(gx*vect, 20*L, label=r"$\mathcal{L}(t)$")
# plt.xlabel(r"$\gamma_x t$", fontsize=17)
# plt.title(r"$h={},\alpha={},j={}$".format(h,alfa,j))
# plt.legend(fontsize=12)
# plt.tight_layout()
# plt.savefig("QSL_LE_Renyi_purestate_j300.pdf")

