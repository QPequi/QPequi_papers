



import numpy as np
from numpy import linalg as LA
from scipy.linalg import expm
# from scipy.special import jv
# from functools import reduce
# from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import time
# from IPython.display import display, Math


import warnings
warnings.filterwarnings('ignore')

plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True



##########################
#
# Dynamics parameters
#
##########################      

j = 250
dim = int(2*j+1)


omega_i = 0
omega_f = np.linspace(0, 2, 50)
w_i = 0
w_f = 0
n_krylov = dim
ti = 0
tf = 100
n_times = 200



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

def mean(v1,A,v2):
    v1_dagger = np.transpose(np.conjugate(v1))
    return np.dot(v1_dagger, np.dot(A, v2))


def time_averaged_complexity(
        omega_i,
        omega_f,
        w_i,
        w_f,
        n_krylov,
        ti,
        tf,
        n_times,
        chi_i = -1,
        chi_f = -1,):
    
    """
        subscript 'i' denotes a pre-quench parameter
        
        subscript 'f' denotes a post-quench parameter
        
        n_krylov is the dimension of the Krylov subspace
        
        ti and tf are the initial and final instant of time of the dynamics
        
        n_times is the number of points of the time evolution
    """
    #Initial hamiltonian
    
    def H_LMG(chi, omega, w):
    # chi = -1
    # Omega = 0
    # w = -1
        return -(chi/2/j)*np.dot(Jz(j),Jz(j)) - omega*Jx(j) - w*Jz(j)
    
    
    #Initial density matrix of the system
    E0, V0 = LA.eig(H_LMG(chi_i, omega_i, w_i))
    index = np.argsort(E0) #Get the index list that would sort E and V
    E0 = E0[index]
    V0 = V0[:,index]
    
    psi0 = V0[:,dim-1]  # Initial state
    
    #######################
    #
    # Parameters of the Krylov
    # basis generation
    #
    #######################
    
    a = np.zeros(n_krylov) #Lanczos coefficients
    b = np.zeros(n_krylov)
    
    K = np.zeros((dim, n_krylov))  #It stores the Krylov states in the columns
    
    """
    #######################
    #
    # Calculating the complexity
    #
    #######################
    """

    
    C = np.zeros(len(omega_f))  # Stores the time averaged complexity
    
    vect = np.linspace(ti, tf, n_times)
    

    for OO in range(1, len(omega_f)):
    
        # chi = -1
        # w = -1
        H = H_LMG(chi_f, omega_f[OO], w_f)
        
        ######################
        #
        # Generating the Krylov basis
        #
        ######################
        
        K[:,0] = psi0  #Lanczos algorithm initiates with the pre-quench ground state
        a[0] = mean(K[:,0], H, K[:,0])
    
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
        
        
        C_t = np.zeros(len(vect))
        
        for tt in range(len(vect)):
            
            # psi_depois = np.dot(chebyshev(H, vect[1]-vect[0], Emax, Emin), psi_antes)
            psit = np.dot(expm(-1j*H*vect[tt]), psi0)
            
            
            for n in range(n_krylov):
                C_t[tt] = C_t[tt] + n*np.abs(np.dot(psit, K[:,n]))**2
            
            
        C[OO] = (1/np.max(vect))*np.trapz(C_t, vect)
        
        # display(Math(r"$\Omega_f = {}$".format(omega_f[OO])))
        print("Omega_f = {}".format(omega_f[OO]))
        
    return(C)


start = time.time()

C = time_averaged_complexity(omega_i,
                             omega_f,
                             w_i,
                             w_f,
                             n_krylov,
                             ti,
                             tf,
                             n_times)


plt.plot(omega_f, C, linewidth = 2)
plt.tick_params(length = 1)
plt.axvline(0.5)

# plt.xlim([0, 2])
# plt.ylim([0, 8.2])


end = time.time()

print("Total time simulation: {} minutes".format((end-start)/60))


