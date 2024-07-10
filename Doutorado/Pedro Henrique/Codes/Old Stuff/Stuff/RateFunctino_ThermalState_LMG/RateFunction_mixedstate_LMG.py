


import numpy as np
from scipy.linalg import expm
#from scipy import linalg as LA
import matplotlib.pyplot as plt
from scipy.signal import find_peaks #it helps to identify t_c (peaks)


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 500
dim = np.int(2*j+1)
gx = 1
gy = 0
vect = np.linspace(0, 30, 300)



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

    

#Initial LMG hamiltonian
h=0
H0 = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))


# #Initial density matrix of the system
# E, V = LA.eig(H0)
# index = np.argsort(E) #Get the index list that would sort E and V
# E = E[index]
# V = V[:,index]


beta = 1

#Initial density matrix: ground state/thermal state
# psi0 = np.dot(V[:,0].reshape(dim,1), V[:,0].reshape(1,dim))
psi0 = expm(-beta*H0)/np.trace(expm(-beta*H0))


#Final LMG hamiltonian
h = 1.2
Hf = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))


######################################################
#
# Rate function for thermal initial state
#
######################################################


R = np.zeros(len(vect)) #rate function

for tt in range(len(vect)):
    U = expm(-1j*vect[tt]*Hf)
    
    LE = np.abs(np.trace(psi0*U))**2
    
    R[tt] = -1/(2*j)*np.log(LE)


#find_peaks command helps to identify where the peaks occur.
# tc, _ = find_peaks(LE, height=0)

plt.plot(gx*vect, R, label=r"$r(t)$")
# plt.plot(gx*vect[tc], LE[tc], 'x')  #mark the peaks throughout r(t)
plt.xlabel("$\gamma_x t$", fontsize=17)
# plt.ylabel(r"$r(t)$", fontsize=17)
plt.legend(fontsize=12, loc=1)
plt.tight_layout()


#np.savetxt("RateFunction_j50_beta1_h2.txt", R)





plt



