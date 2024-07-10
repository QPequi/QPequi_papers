
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
# from scipy import integrate


# plt.rcParams["font.size"] = 15
# plt.rcParams["text.usetex"] = True


##########################
#
# Dynamics parameters
#
##########################      

j = 25
dim = int(2*j+1)


omega_i = 0.6
omega_f = np.linspace(0,1,30)
w_i = 0
w_f = [0]
ti = 0
tf = 150
n_times = 300

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

def mean(v1,A,v2):
    v1_dagger = np.transpose(np.conjugate(v1))
    return np.dot(v1_dagger, np.dot(A, v2))


def mag(omega_i,
        omega_f,
        w_i,
        w_f,
        ti,
        tf,
        n_times,
        chi_i = 1,
        chi_f = 1):
    
    """
        subscript 'i' denotes a pre-quench parameter
        
        subscript 'f' denotes a post-quench parameter
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
    
    M_t = np.zeros(n_times)
    vect = np.linspace(ti,tf,n_times)
    M_av = np.zeros(len(omega_f))
    
    for oo in range(len(omega_f)):
        H = H_LMG(chi_f,omega_f[oo],w_f)
        print(omega_f[oo])
        
        for tt in range(n_times):
            
            psit = np.dot(expm(-1j*H*vect[tt]), psi0)
            
            M_t[tt] = np.real(mean(psit, Jz(j), psit))
    
        M_av[oo] = (1/np.max(vect))*np.trapz(M_t, vect)
        
    return(M_av)


M_av = mag(omega_i,
            omega_f,
            w_i,
            w_f,
            ti,
            tf,
            n_times)

# np.savetxt("TimeAveragedMagnetization_vs_h_N{}.txt".format(2*j), M_av)

plt.plot(omega_f, M_av/j)
# plt.legend(fontsize=10, loc=1)  






