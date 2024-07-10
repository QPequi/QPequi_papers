

#########################
#
# Inverse participation ratio
# in the Krylov basis
#
#########################



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



##########################
#
# Dynamics parameters
#
##########################      



def IPR_Krylov(j,
               omega_i,
               omega_f,
               w_i,
               w_f,
               n_krylov,
               ti,
               tf,
               n_times,
               basis,
               chi_i = 1,
               chi_f = 1):
            
            
            dim = int(2*j+1)
            
            
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
            # 0 -> all spins pointing up
            # 1 -> all spins pointing down
            
            #Pos-quench hamiltonian
            H = H_LMG(chi_f, omega_f, w_f)
            
            
            ######################
            #
            # Generating the Krylov basis
            #
            ######################
            
            # if basis =='krylov':
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
            # Calculating IPR
            #
            ######################
            
            vect = np.linspace(ti, tf, n_times)
            IPR_t = np.zeros(len(vect))
            
            for tt in range(len(vect)):
                
                print(tt)
                U = expm(-1j*H*vect[tt])
                
                for k in range(n_krylov):
                    
                    if basis == 'energy':
                        IPR_t[tt] = IPR_t[tt] + np.abs(expval(V0[:,k], U, psi0))**4
                    
                    elif basis == 'krylov':
                        IPR_t[tt] = IPR_t[tt] + np.abs(expval(K[:,k], U, psi0))**4
            
            
            return(IPR_t)
 
        
j = 100
#omega_i = 0.3
omega_f = 0.8
w_i = 0
w_f = [0]
n_krylov = int(2*j+1)
ti = 0
tf = 50
n_times = 350
vect = np.linspace(ti,tf,n_times)


for basis in ['krylov', 'energy']: #distinguish the initial degenerated states
    i = 0
    for omega_i in [0.4]:
        
        
        IPR_t = IPR_Krylov(j,
                        omega_i,
                        omega_f,
                        w_i,
                        w_f,
                        n_krylov,
                        ti,
                        tf,
                        n_times,
                        basis)
        
        # np.savetxt("Krylov_IPR_LMG_h{}_N{}_d{}.txt".format(omega_f[0], 2*j,d), IPR_t)
        
        
        if omega_f <= 1/2:
            c = ['#003366', '#0080FF', '#99CCFF']  #blue tones
        elif omega_f > 1/2:
            c = ['#660000', '#FF0000', '#FF9999']  #red tones
        
        
        # IPR_t = np.loadtxt("Krylov_IPR_LMG_h0.3_N200_d{}.txt".format(d))
        
        # Plotting
        if basis=='energy':
            plt.plot(vect, IPR_t, '--', color=c[i+1])
            plt.xlabel(r"$t$", fontsize=20)
            # plt.ylabel(r"$IPR^+(t)/j, C^-(t)$", fontsize=20)
            # plt.legend(frameon=False, fontsize=16, 
                        # loc=(0.6,0.72))
                        # loc=1)
        
        elif basis=='krylov':
            plt.plot(vect, IPR_t, '-.', color=c[i])
            plt.xlabel(r"$t$", fontsize=20)
            plt.ylabel(r"$IPR(t)$", fontsize=20)
        
        i = i+1
    
    plt.text(38.5, 0.9, r"$h_0={}$".format(omega_i), fontsize=18)
    plt.legend(frameon=False, fontsize=16, 
                    # loc=(0.6,0.72))
                    loc=1)
    plt.tight_layout()
            
            
        
plt.savefig("KrylovPR_vs_t_LMG_h0{}_hf{}_N{}_1.pdf".format(omega_i,omega_f,2*j))

