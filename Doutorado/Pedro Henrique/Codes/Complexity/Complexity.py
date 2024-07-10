

#######################################
#
# This code calculates the evolution of
# the quantum state spread complexity under
# a quench dynamics in the LMG model
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


    
    
def complexity(j,
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
    
    psi0 = V0[:,1]  # Initial state
    # d=0 -> collective spin pointing along +z
    # d=1 -> collective spin pointing along -z
    
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
    
    K, R = LA.qr(K)
        
    ######################
    #
    # Calculating the complexity
    #
    ######################
    
    vect = np.linspace(ti, tf, n_times)
    C_t = np.zeros(len(vect))
    
    for tt in range(len(vect)):
        
        print(tt)
        # psi_depois = np.dot(chebyshev(H, vect[1]-vect[0], Emax, Emin), psi_antes)
        psit = np.dot(expm(-1j*H*vect[tt]), psi0)
                
                
        for n in range(n_krylov):
            C_t[tt] = C_t[tt] + n*np.abs(np.dot(psit, K[:,n]))**2
            
    return(C_t, K)
    

##########################
#
# Dynamics parameters
#
##########################      


w_i = 0
w_f = [0]
ti = 0
tf = 150
n_times = 350
vect = np.linspace(ti,tf,n_times)


C_t = []

for j in [100]:
    for omega_i, omega_f in [(0.1, 0.4), (0.2, 0.4), (0.1, 0.8), (0.2, 0.8)]:
        n_krylov = int(2*j+1)
        
        C_t.append(complexity(j,
                       omega_i,
                       omega_f,
                       w_i,
                       w_f,
                       n_krylov,
                       ti,
                       tf,
                       n_times)[0])
        
        

if omega_f <= 1/2:
    c = ['#003366', '#0080FF', '#99CCFF']  #blue tones
elif omega_f > 1/2:
    c = ['#660000', '#FF0000', '#FF9999']  #red tones


# Plotting


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

ax1.plot(vect, C_t[0]/j)
ax1.text(130, 0, r'$(c)$')
ax1.text(88, 0.19, r'$h_0={}$'.format(0.1), fontsize=18)
ax1.text(88, 0.15, r'$h_f={}$'.format(0.4), fontsize=18)
ax1.set_ylabel(r'$C_\mathcal{K}(t)/j$', fontsize=20)

ax2.plot(vect, C_t[1]/j)
ax2.text(130, 0, r'$(d)$')
ax2.text(88, 0.082, r'$h_0={}$'.format(0.2), fontsize=18)
ax2.text(88, 0.062, r'$h_f={}$'.format(0.4), fontsize=18)
# ax2.set_ylim([0, 0.25])

ax3.plot(vect, C_t[2]/j, color='darkred')
ax3.plot(vect, np.ones(len(vect)), color='black', linestyle='--', linewidth=1)
ax3.text(130, 0.1, r'$(e)$')
ax3.text(92, 1.65, r'$h_0={}$'.format(0.1), fontsize=18)
ax3.text(92, 1.3, r'$h_f={}$'.format(0.8), fontsize=18)
ax3.set_xlabel(r'$t$', fontsize=20)
ax3.set_ylabel(r'$C_\mathcal{K}(t)/j$', fontsize=20)

ax4.plot(vect, C_t[3]/j, color='darkred')
ax4.plot(vect, np.ones(len(vect)), color='black', linestyle='--', linewidth=1)
ax4.text(130, 0.1, r'$(f)$')
ax4.text(92, 1.65, r'$h_0={}$'.format(0.2), fontsize=18)
ax4.text(92, 1.3, r'$h_f={}$'.format(0.8), fontsize=18)
ax4.set_xlabel(r'$t$', fontsize=20)
ax4.set_ylim([-0.05, 2.01])
plt.tight_layout(pad=0.3)

# for ax in fig.get_axes():
#     ax.label_outer()
    

plt.savefig("Complexity_vs_t_LMG_N{}_2.pdf".format(2*j))

