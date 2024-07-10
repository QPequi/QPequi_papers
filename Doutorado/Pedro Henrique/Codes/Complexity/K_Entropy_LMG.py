


#######################################
#
# This code calculates the evolution of
# the K-entropy for the quenched LMG model
#
#######################################

import numpy as np
from numpy import linalg as LA
from scipy.linalg import expm
# from functools import reduce
# from scipy.signal import find_peaks
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

plt.rcParams["font.size"] = 22
plt.rcParams["text.usetex"] = True

   


        
def K_Entropy(j,
            omega_i,
            omega_f,
            w_i,
            w_f,
            n_krylov,
            ti,
            tf,
            n_times,
            chi_i = 1,
            chi_f = 1):
        
    
        dim = int(2*j+1)
        
        #Defining the Jx, Jy and Jz in z basis and other useful functions
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
        
        def expval(v1,A,v2):
            v1_dagger = np.transpose(np.conjugate(v1))
            return np.dot(v1_dagger, np.dot(A, v2))
            
        # LMG hamiltonian
        def H_LMG(chi, omega, w):
            return -(chi/(2*j))*np.dot(Jz(j),Jz(j)) - omega*Jx(j) - w*Jz(j)
        
    
        #Pre-quench eigenvalues and eigenvectors
        E0, V0 = LA.eig(H_LMG(chi_i, omega_i, w_i))
        index = np.argsort(E0) #Get the index list that sorts E and V
        E0 = E0[index]
        V0 = V0[:,index]
        
        # d=1
        psi0 = V0[:,1]  # Initial state
        # d=0 -> collective spin pointing up
        # d=1 -> collective spin pointing down
        
        #Pos-quench hamiltonian
        H = H_LMG(chi_f, omega_f, w_f)
        
        
        ######################
        #
        # Generating the Krylov basis
        #
        ######################
        
        
        a = np.zeros(n_krylov) #Lanczos coefficients
        b = np.zeros(n_krylov)
        
        K = np.zeros((dim, n_krylov))  #It stores the Krylov states in the columns
        
        K[:,0] = psi0  #Lanczos algorithm initiates with the pre-quench ground state
        a[0] = expval(K[:,0], H, K[:,0])
    
        K[:,1] = np.dot(H, K[:,0]) - a[0]*K[:,0]
        b[1] = LA.norm(K[:,1])  #normalizing K_1
        K[:,1] = K[:,1]/b[1]
    
        for s in range(2, n_krylov):  #Generating the next Krylov states
            a[s-1] = np.dot(K[:,s-1], np.dot(H, K[:,s-1]))
            K[:,s] = np.dot(H, K[:,s-1]) - a[s-1]*K[:,s-1] - b[s-1]*K[:,s-2] #Lanczos algorithm
    
            b[s] = LA.norm(K[:,s]) #Normalization
            K[:,s] = K[:,s]/b[s]
            
        
        #re-orthogonalization applying modified Gram-Schmidt explicitly
        
        # for n in range(n_krylov):
            
        #     b[n] = np.linalg.norm(K[:, n])
        #     K[:, n] = K[:, n] / b[n]
            
        #     for m in range(n+1, n_krylov):
        #         K[:,n] = K[:,n] - np.dot(K[:, m], K[:, n]) * K[:, m]
        
        
        K, R = np.linalg.qr(K)
        
        ######################
        #
        # Calculating the complexity
        #
        ######################
        
        vect = np.linspace(ti, tf, n_times)
        H_Energy = np.zeros(len(vect))
        H_Krylov = np.zeros(len(vect))
        
        for tt in range(len(vect)):
            
            print(tt)
            # psi_depois = np.dot(chebyshev(H, vect[1]-vect[0], Emax, Emin), psi_antes)
            psit = np.dot(expm(-1j*H*vect[tt]), psi0)
                    
                    
            for n in range(n_krylov):
                pK_n = np.abs(np.dot(psit, K[:,n]))**2
                pE_n = np.abs(np.dot(psit, V0[:,n]))**2
                
                if pK_n != 0:
                    
                    H_Krylov[tt] = H_Krylov[tt] - pK_n*np.log(pK_n)
                
                if pE_n != 0:
                    H_Energy[tt] = H_Energy[tt] - pE_n*np.log(pE_n)
                
        return(H_Energy, H_Krylov)
  

j = 100
omega_i = 0.1
omega_f = 0.4
w_i = 0
w_f = 0
n_krylov = int(2*j+1)
ti = 0
tf = 250
n_times = 350      
vect = np.linspace(ti,tf,n_times)


H_Energy = []
H_Krylov = []

for j in [100]:
    for omega_i, omega_f in [(0.1, 0.4), (0.2, 0.4), (0.1, 0.8), (0.2, 0.8)]:
        n_krylov = int(2*j+1)

        H = K_Entropy(j,
                     omega_i,
                     omega_f,
                     w_i,
                     w_f,
                     n_krylov,
                     ti,
                     tf,
                     n_times)
                
        H_Energy.append(H[0])
        H_Krylov.append(H[1])


#################
#
# Plot
#
#################


if omega_f <= 1/2:
    c = ['#003366', '#0080FF', '#99CCFF']  #blue tones
elif omega_f > 1/2:
    c = ['#660000', '#FF0000', '#FF9999']  #red tones

i=0

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

ax1.plot(vect, H_Energy[0]/j, color='#0080FF')
ax1.plot(vect, H_Krylov[0]/j, color='#003366', linestyle='--')


ax1.text(217, 0.0018, r'$(a)$')
ax1.text(160, 0.018, r'$h_0={}$'.format(0.1), fontsize=15)
ax1.text(160, 0.012, r'$h_f={}$'.format(0.4), fontsize=15)
ax1.set_ylabel(r'$\mathcal{E}(t)/j$', fontsize=20)
# ax1.set_ylim([-0.002, 0.035])

##############
##############

ax2.plot(vect, H_Energy[1]/j, color='#0080FF')
ax2.plot(vect, H_Krylov[1]/j, color='#003366', linestyle='--')

ax2.text(220, 0.0018, r'$(b)$')
ax2.text(160, 0.015, r'$h_0={}$'.format(0.2), fontsize=15)
ax2.text(160, 0.01, r'$h_f={}$'.format(0.4), fontsize=15)
# ax2.set_ylim([-0.002, 0.035])


##############
##############

ax3.plot(vect, H_Energy[2]/j, color='#FF0000')
ax3.plot(vect, H_Krylov[2]/j, color='#660000', linestyle='--')

# ax3.plot(vect, np.ones(len(vect)), color='black', linestyle='--', linewidth=1)
ax3.text(220, 0.0018, r'$(c)$')
ax3.text(160, 0.022, r'$h_0={}$'.format(0.3), fontsize=15)
ax3.text(160, 0.014, r'$h_f={}$'.format(0.8), fontsize=15)
ax3.set_xlabel(r'$t$', fontsize=20)
ax3.set_ylabel(r'$\mathcal{E}(t)/j$', fontsize=20)


#############
#############


ax4.plot(vect, H_Energy[3]/j, color='#FF0000')
ax4.plot(vect, H_Krylov[3]/j, color='#660000', linestyle='--')
# ax4.plot(vect, np.ones(len(vect)), color='black', linestyle='--', linewidth=1)
ax4.text(217, 0.0018, r'$(d)$')
ax4.text(160, 0.022, r'$h_0={}$'.format(0.4), fontsize=15)
ax4.text(160, 0.014, r'$h_f={}$'.format(0.8), fontsize=15)
ax4.set_xlabel(r'$t$', fontsize=20)
# ax4.set_ylim([-0.05, 2.01])


plt.tight_layout(pad=0.3)


plt.savefig("KEntropy_vs_t_LMG_N{}.eps".format(2*j))
# plt.savefig("ab.pdf")
