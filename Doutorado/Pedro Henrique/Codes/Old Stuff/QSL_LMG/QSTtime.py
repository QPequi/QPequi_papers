



import numpy as np
from numpy import linalg as LA
from scipy.linalg import fractional_matrix_power
from scipy.linalg import expm
# from scipy.linalg import logm
from scipy.signal import find_peaks
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 100
dim = np.int(2*j+1)
gx = 1
gy = 0
vect = np.linspace(0, 10, 100)



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

#Commutator 
def Comm(A,B):
    return np.dot(A,B)-np.dot(B,A)


#Schatten 2-norm
def norm2(A):
    return np.sqrt(np.trace(np.dot(np.transpose(np.conjugate(A)), A)))


#Finding minimuns
def min(A):
    M = np.zeros(len(vect))
    m = []
    M = np.r_[True, A[1:] < A[:-1]] & np.r_[A[:-1] < A[1:], True]
    
    for i in range(len(M)):
        if M[i]==True:
            m.append(i)
    return m


#Initial LMG hamiltonian
h=0
H0 = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))


#Initial density matrix of the system
E, V = LA.eig(H0)
index = np.argsort(E) #Get the index list that would sort E and V
E = E[index] 
V = V[:,index]

# beta = 1

#Initial density matrix: ground state/thermal state
rho0 = np.dot(V[:,0].reshape(dim, 1), V[:,0].reshape(1, dim))
# rho0 = expm(-beta*H0)/np.trace(expm(-beta*H0))



E0, V0 = LA.eig(rho0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index] 
V0 = V0[:,index]


Tc = np.zeros(40) # first maximum time of QSL
i=0
for h in [0.85]:
        
    #Final LMG hamiltonian
    Hf = -h*Jz(j) - (gx/2/j)*np.dot(Jx(j),Jx(j)) - (gy/2/j)*np.dot(Jy(j),Jy(j))
    
    
    #Final equilibrium state 
    # rhof_eq = expm(-beta*Hf)/np.trace(expm(-beta*Hf))
    
    # rhof_eig, rhof_v = LA.eig(rhof_eq)
    
    ####################################################
    #
    # Relative entropy 
    #
    ####################################################
    
    
    
    RRE = np.zeros(len(vect)) #Rényi relative entropy
    # TRE = np.zeros(len(vect)) #Tsallis relative entropy
    
    tau = np.zeros(len(vect)) #QSL time
    
    
    alpha = 0.6
        
    G = np.abs(1+(1-alpha)*np.log(np.min(E0)))**(-1)*norm2(fractional_matrix_power(rho0, 1-alpha))*norm2(Comm(Hf, fractional_matrix_power(rho0, alpha)))
    # G = norm2(fractional_matrix_power(rho0, 1-alfa))*norm2(Comm(Hf, fractional_matrix_power(rho0, alfa)))
        
    
    for tt in range(len(vect)):
        U = expm(-1j*Hf*vect[tt])  #time evolution operator
        Ut = np.transpose(np.conjugate(U))
        
        rhot = np.dot(U, np.dot(rho0, Ut))  #quenched density matrix at time tt
        
        RRE = 1/(alpha-1)*np.log(np.trace(np.dot(fractional_matrix_power(rhot, alpha), fractional_matrix_power(rho0, 1-alpha))))
        # TRE = 1/(1-alfa)*(1-np.trace(np.dot(fractional_matrix_power(rhot, alfa), fractional_matrix_power(rho0, 1-alfa))))
        
    
        tau[tt] = np.abs(1-alpha)*RRE/G
    
    
    #find_peaks command helps to find peaks
    tc, _ = find_peaks(tau, height=0)
    
    Tc[i] = np.int(tc[0])
    i = i+1

Tc = [np.int(x) for x in Tc]

# plt.plot(np.arange(0.01, 2, 0.05), Tc, 'o-', label=r"$\alpha={}, j={}$".format(alpha,j))
# plt.xlabel(r"$h_f$", fontsize=17)
# plt.ylabel(r"$T_0$", fontsize=17)
# plt.legend(fontsize=12)
# plt.xlim([0,2.001])
# plt.grid(True)
# plt.tight_layout()
# plt.savefig("FirstMaxQSL_vs_hf_j100_alpha0.6.pdf")


plt.plot(gx*vect, tau, 'o-', label=r"$h={},\alpha={},j={}$".format(h,alpha,j))
plt.plot(gx*vect[tc], tau[tc], '*')
# plt.plot(vect, -150/(2*j)*np.log(L), label='$h={}, \gamma_x={}, j={}$'.format(h,gx,j))
# plt.plot(vect, L, label='$h={}, \gamma_x={}, j={}$'.format(h,gx,j))
# plt.xlabel("$\gamma_x t$", fontsize=17)
# # plt.xlim([0,30])
# plt.ylabel(r"$\tau^{R}_{\alpha}(\rho(t)||\rho_0)$", fontsize=17)
# # plt.title(r"$h={},j={}$".format(h,j))
# plt.legend(fontsize=12)
# plt.tight_layout()


# plt.savefig("QSL_Renyi_purestate_j300_h0.5.pdf")



