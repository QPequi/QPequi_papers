




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
# import time

import warnings
warnings.filterwarnings('ignore')

plt.rcParams["font.size"] = 22
plt.rcParams["text.usetex"] = True



##########################
#
# Dynamics parameters
#
##########################      


k = 0
for j in [100]:
    for omega_f in [0.3, 0.5, 0.7]:
        dim = int(2*j+1)
        
        
        omega_i = 0
        # omega_f = [0.8]
        w_i = 0
        w_f = [0]
        n_krylov = dim+1
        ti = 0
        tf = 250
        n_times = 350
        
        
        
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
            
        
        
        def vN_Entropy(
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
            
            
            # LMG hamiltonian
            def H_LMG(chi, omega, w):
                return -(chi/(2*j))*np.dot(Jz(j),Jz(j)) - omega*Jx(j) - w*Jz(j)
            
        
            #Pre-quench eigenvalues and eigenvectors
            E0, V0 = LA.eig(H_LMG(chi_i, omega_i, w_i))
            index = np.argsort(E0) #Get the index list that sorts E and V
            E0 = E0[index]
            V0 = V0[:,index]
            
            # d=0
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
                
            
            ######################
            #
            # Calculating the complexity
            #
            ######################
            
            # vect = np.linspace(ti, tf, n_times)
            # H_Energy = np.zeros(len(vect))
            # H_Krylov = np.zeros(len(vect))
            
            # for tt in range(len(vect)):
                
            #     print(tt)
            #     # psi_depois = np.dot(chebyshev(H, vect[1]-vect[0], Emax, Emin), psi_antes)
            #     psit = np.dot(expm(-1j*H*vect[tt]), psi0)
                        
                        
            #     for n in range(n_krylov):
            #         pK_n = np.abs(np.dot(psit, K[:,n]))**2
            #         pE_n = np.abs(np.dot(psit, V0[:,n]))**2
                    
            #         if pK_n != 0:
                        
            #             H_Krylov[tt] = H_Krylov[tt] - pK_n*np.log(pK_n)
                    
            #         if pE_n != 0:
            #             H_Energy[tt] = H_Energy[tt] - pE_n*np.log(pE_n)
                    
            return(V0, K, E0, b)
        
        vect = np.linspace(ti,tf,n_times)
        
        V0, K, E0, b = vN_Entropy(omega_i,
                       omega_f,
                       w_i,
                       w_f,
                       n_krylov,
                       ti,
                       tf,
                       n_times)
        
        
        if omega_f <= 1/2:
            c = ['#003366', '#0080FF', '#99CCFF']  #blue tones
        elif omega_f > 1/2:
            c = ['#660000', '#FF0000', '#FF9999']  #red tones
        
        plt.plot(b, 'o', color=c[k], label=r"$h_f={}$".format(omega_f))
        plt.ylabel(r"$b_{m_z}$", fontsize=22)
        plt.xlabel(r"$m_z$", fontsize=22)
        # plt.legend(frameon=False, fontsize=10, loc=1)
        # plt.tight_layout()
        
        mz = np.arange(0, 2*j+2, 1)
        b_theor = np.zeros(len(mz))
        
        for i in range(len(mz)):
            b_theor[i] = (omega_f/2)*np.sqrt(mz[i]*(2*j-mz[i]+1))
            
        
        plt.plot(b_theor, color='yellow')
        plt.legend(frameon=False, fontsize=16, loc=(0.3, 0.01))
        plt.tight_layout()
        k += 1

plt.savefig("LanczosCoeff_N{}_+.eps".format(2*j))
        
        # if omega_f <= 1/2:
        #     c = ['#003366', '#0080FF', '#99CCFF']  #blue tones
        # elif omega_f > 1/2:
        #     c = ['#660000', '#FF0000', '#FF9999']  #red tones
        
        
        
        # plt.plot(vect, H_Krylov/j, '--', color=c[i], label=r"$Krylov\ basis$")
        # plt.xlabel(r"$t$", fontsize=20)
        # # plt.ylabel(r"$H_E(t)/j$", fontsize=20)
        # plt.legend(frameon=False, fontsize=16, 
        #            #loc=(0.6,0.72))
        #           loc=1)
        
        
        # plt.plot(vect, H_Energy/j, '-.', color=c[i], label=r"$Energy\ basis$")
        # plt.xlabel(r"$t$", fontsize=20)
        # plt.ylabel(r"$H(t)/j$", fontsize=20)
        # plt.legend(frameon=False, fontsize=16)
        
        # # plt.annotate('$N={}$'.format(2*j), xytext=(1,1), fontsize=16)
        
        # plt.tight_layout()
        # i = i+1

# plt.savefig("Entropy_EandK_vs_t_LMG_h{}_N{}.pdf".format(omega_f,2*j))

