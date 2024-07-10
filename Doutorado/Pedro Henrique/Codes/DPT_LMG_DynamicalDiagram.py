

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
# from scipy import integrate
import time


import warnings
warnings.filterwarnings('ignore')

plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


##########################
#
# Dynamics parameters
#
##########################      



j = 50
dim = int(2*j+1)
n_pixels = 50

for omega_i in [0, 0.2, 0.4]:
    # omega_i = 0
    print(omega_i)
    omega_f = np.linspace(0, 1, n_pixels)
    w_i = 0
    w_f = [0]
    n_krylov = dim
    ti = 0
    tf = 150
    n_times = 350
    
    
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
    
    def mean(v1,A,v2):
        v1_dagger = np.transpose(np.conjugate(v1))
        return np.dot(v1_dagger, np.dot(A, v2))
    
    def diagram_magnetization(
            omega_i,
            omega_f,
            w_i,
            w_f,
            ti,
            tf,
            n_times,
            chi_i = 1,
            chi_f = 1):
        
        """
            subscript 'i' identifies a pre-quench parameter
            
            subscript 'f' identifies a post-quench parameter
            [w_f is expected to be a list or an array]
            
            ti and tf are the initial and final instant of time of the dynamics
            
            n_times is the number of points of the time evolution
        """
        #Initial hamiltonian
        
        def H_LMG(chi, omega, w):
    
            return -(chi/(2*j))*np.dot(Jz(j),Jz(j)) - omega*Jx(j) - w*Jz(j)
        
    
        #Initial density matrix of the system
        E0, V0 = LA.eig(H_LMG(chi_i, omega_i, w_i))
        index = np.argsort(E0) #Get the index list that would sort E and V
        E0 = E0[index]
        V0 = V0[:,index]
        
        psi0 = V0[:,1]  # Initial state 
        
        
        #######################
        #
        # Caculating the magnetization
        #
        #######################
        
        vect = np.linspace(ti, tf, n_times)
        
        M = np.zeros((len(omega_f) , len(w_f)))  # Stores the time averaged complexity
        M_t = np.zeros(len(vect))
        vect = np.linspace(ti, tf, n_times)
        
    
        for OO in range(len(omega_f)-1, -1, -1): #run backwards through the values of omega_f
            for ww in range(len(w_f)):
        
                H = H_LMG(chi_f, omega_f[OO], w_f[ww])
                
                    
                ######################
                #
                # Calculating the magnetization
                #
                ######################
                
                
                for tt in range(len(vect)):
                    
                    # psi_depois = np.dot(chebyshev(H, vect[1]-vect[0], Emax, Emin), psi_antes)
                    psit = np.dot(expm(-1j*H*vect[tt]), psi0)
                    
                    M_t[tt] = np.dot(np.conjugate(psit), np.dot(Jz(j), psit))
                    
                    
                M[len(omega_f)-1-OO,ww] = (1/np.max(vect))*np.trapz(M_t, vect)
                print("(omega_f, w_f) = ({},{})".format(omega_f[OO], w_f[ww]))
                
        return(M)
    
    
    start = time.time()
    
    M = diagram_magnetization(omega_i,
                            omega_f,
                            w_i,
                            w_f,
                            ti,
                            tf,
                            n_times)
    
    
    
    if len(w_f) > 1:
        plt.imshow(M/j, cmap='rainbow',
                    extent =[w_f.min(), w_f.max(), omega_f.min(), omega_f.max()], 
                    aspect='auto')
        plt.colorbar(orientation="vertical")
        plt.xlabel(r"$\omega$", fontsize=17)
        plt.ylabel(r"$\Omega$", fontsize=17)
        plt.title(r"$2\overline{S_z}/N$", loc='right')
        plt.tick_params(length = 1)
        plt.tight_layout()
    #     # np.savetxt("2DDiagram_MagnetizationSz_j50.txt", M)
    #     # plt.savefig("2D_diagram_TimeAveMagnetization_N200.pdf")
    
    else:
        M = M[:,0]
        plt.plot(omega_f, M[::-1]/j, label=r"$N={}$".format(2*j))
        plt.xlabel(r"$\Omega$", fontsize=17)
        plt.ylabel(r"$\overline{S_z}$", fontsize=17)
        plt.legend(frameon=False, fontsize=12)
        plt.tight_layout()

# plt.savefig("TimeAveragedComplexity_w0.pdf")



# end = time.time()









