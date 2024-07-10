



import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
# from scipy import integrate
import time


import warnings
warnings.filterwarnings('ignore')

plt.rcParams["font.size"] = 22
plt.rcParams["text.usetex"] = True


##########################
#
# Dynamics parameters
#
##########################      

i=0
for j in [100]:
    dim = int(2*j+1)
    
    omega_i = 0.1
    omega_f = 0.4
    w_i = 0
    w_f = [0]
    ti = 0
    tf = 250
    n_times = 350
    
    vect = np.linspace(ti, tf, n_times)
    
    
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
    
    
    def magnetization(
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
        
        
        psi0 = V0[:,1]
        #psi0 = V0[:,0]  # Initial state 
        
            
        ######################
        #
        # Calculating the magnetization
        #
        ######################
        
        vect = np.linspace(ti, tf, n_times)
        
        M_t = np.zeros(len(vect))
        vect = np.linspace(ti, tf, n_times)
        
        
        H = H_LMG(chi_f, omega_f, w_f)
        
        for tt in range(len(vect)):
            
            print(tt)
            # psi_depois = np.dot(chebyshev(H, vect[1]-vect[0], Emax, Emin), psi_antes)
            psit = np.dot(expm(-1j*H*vect[tt]), psi0)
            
            M_t[tt] = np.dot(np.conjugate(psit), np.dot(Jz(j), psit))
            
                
        return(M_t, psit, psi0, E0, V0)
    
    
    start = time.time()
    
    M_t, psit,psi0, E0, V0 = magnetization(omega_i,
                        omega_f,
                        w_i,
                        w_f,
                        ti,
                        tf,
                        n_times)
    
    if omega_f <= 1/2:
        c = ['#003366', '#0080FF', '#99CCFF']  #blue tones
    elif omega_f > 1/2:
        c = ['#660000', '#FF0000', '#FF9999']  #red tones
        
    plt.plot(vect, M_t/j, color=c[i], label=r"$N={}$".format(2*j))
    plt.xlabel(r"$t$", fontsize=20)
    plt.ylabel(r"$S_z(t)/j$", fontsize=20)
    plt.legend(frameon=False, fontsize=16)
               # loc=(0.55, 0))
               
    # plt.plot(vect, M_t/j, '-.', color=c[i], label=r"$C^-(t)/j$")
    # plt.xlabel(r"$t$", fontsize=20)
    # # plt.ylabel(r"$C^-(t)/j, C^-(t)$", fontsize=20)
    # plt.legend(frameon=False, fontsize=16, 
    #             #loc=(0.6,0.72))
    #             loc=1)
    
    plt.tight_layout()
    # plt.savefig("TimeAveragedComplexity_w0_N50.pdf")

    i += 1
    end = time.time()



# plt.savefig("Magnetization_vs_t_h{}.pdf".format(omega_f[0]))





