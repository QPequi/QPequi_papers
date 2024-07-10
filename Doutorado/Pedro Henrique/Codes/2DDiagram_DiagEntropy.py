




import numpy as np
from numpy import linalg as LA
from scipy.linalg import expm
# from functools import reduce
# from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import time

import warnings
warnings.filterwarnings('ignore')

# plt.rcParams["font.size"] = 15
# plt.rcParams["text.usetex"] = True


##########################
#
# Dynamics parameters
#
##########################      

j = 25
dim = int(2*j+1)


omega_i = 0
omega_f = np.linspace(0, 1, 51)
w_i = 0
w_f = np.linspace(-1,1,51)
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


def diagram_diagentropy(
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
    
    """
        subscript 'i' denotes a pre-quench parameter
        
        subscript 'f' denotes a post-quench parameter
        [w_f is expected to be a list or an array]
        
        n_krylov is the dimension of the Krylov subspace
        
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
    
    # psi0 = V0[:,1]  # Initial state 
    
    
    # #######################
    # #
    # # Parameters of the Krylov
    # # basis generation
    # #
    # #######################

    # a = np.zeros(n_krylov) #Lanczos coefficients
    # b = np.zeros(n_krylov)
    
    # K = np.zeros((dim, n_krylov))  #It stores the Krylov states in the columns
    
    
    #######################
    #
    # Calculating the Diagonal entropy
    #
    #######################

    
    S = np.zeros((len(omega_f) -1, len(w_f)))  # Stores the time averaged complexity
    
    vect = np.linspace(ti, tf, n_times)
    

    for OO in range(len(omega_f)-1, 0, -1): #run backwards through the values of omega_f
        for ww in range(len(w_f)):
    
            H = H_LMG(chi_f, omega_f[OO], w_f[ww])
            
            # ######################
            # #
            # # Generating the Krylov basis
            # #
            # ######################
            
            # K[:,0] = psi0  #Lanczos algorithm initiates with the pre-quench ground state
            # a[0] = mean(K[:,0], H, K[:,0])
        
            # K[:,1] = np.dot(H, K[:,0]) - a[0]*K[:,0]
            # b[1] = LA.norm(K[:,1])  #normalizing K_1
            # K[:,1] = K[:,1]/b[1]
        
            # for s in range(2, n_krylov):  #Generating the next Krylov states
            #     a[s-1] = np.dot(K[:,s-1], np.dot(H, K[:,s-1]))
            #     K[:,s] = np.dot(H, K[:,s-1]) - a[s-1]*K[:,s-1] - b[s-1]*K[:,s-2] #Lanczos algorithm
        
            #     b[s] = LA.norm(K[:,s]) #Normalization
            #     K[:,s] = K[:,s]/b[s]
                
            ######################
            #
            # Calculating the complexity
            #
            ######################
            
            
            S_t = np.zeros(len(vect))
            
            for tt in range(len(vect)):
                U = expm(-1j*H*vect[tt])
                
                for k in range(int(2)):
                    Ck = np.abs(np.dot(np.conjugate(V0[:,k]), np.dot(U, V0[:,dim-1])))**2
                    
                    if Ck < 9e-5:
                        S_t[tt] = S_t[tt] + 0
                        
                    else:
                        S_t[tt] = S_t[tt] - Ck*np.log(Ck)
                
                
            S[len(omega_f)-1-OO,ww] = (1/np.max(vect))*np.trapz(S_t, vect)
            
            # display(Math(r"$\Omega_f = {}$".format(omega_f[OO])))
            print("(omega_f, w_f) = ({},{})".format(omega_f[OO], w_f[ww]))
            
    return(S)

start = time.time()

S = diagram_diagentropy(omega_i,
                        omega_f,
                        w_i,
                        w_f,
                        n_krylov,
                        ti,
                        tf,
                        n_times)


if len(w_f) > 1:
    plt.imshow(S/j, cmap='rainbow',
                extent =[w_f.min(), w_f.max(), omega_f.min(), omega_f.max()], 
                aspect='auto')
    plt.colorbar(orientation="vertical")
    plt.xlabel(r"$\omega$", fontsize=17)
    plt.ylabel(r"$\Omega$", fontsize=17)
    plt.title(r"$2\overline{C}/N$", loc='right')
    plt.tick_params(length = 1)
    plt.tight_layout()
    # plt.savefig("2D_diagram_TimeAveragedComplexity_N50.pdf")

else:
    S = S[:,0]
    plt.plot(omega_f[1:], S[::-1]/j, '-o', label=r"$\omega = 0$")
    plt.xlabel(r"$\Omega$", fontsize=17)
    plt.ylabel(r"$\overline{C}$", fontsize=17)
    plt.legend(fontsize=10)
    plt.tight_layout()
    # plt.savefig("TimeAveragedComplexity_w0_N50.pdf")


end = time.time()

# np.savetxt("2DDiagram_Complexity_N1000.pdf", C)

# print("Total time simulation: {} minutes".format((end-start)/60))






