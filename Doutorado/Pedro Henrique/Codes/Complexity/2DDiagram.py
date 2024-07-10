



import numpy as np
from numpy import linalg as LA
from scipy.linalg import expm
from scipy.integrate import simps
# from functools import reduce
# from scipy.signal import find_peaks
import matplotlib.pyplot as plt
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


        
        
def diagram_complexity(j,
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
    #Initial hamiltonian
    
    def H_LMG(chi, omega, w):

        return -(chi/(2*j))*np.dot(Jz(j),Jz(j)) - omega*Jx(j) - w*Jz(j)
    

    #Initial density matrix of the system
    E0, V0 = LA.eig(H_LMG(chi_i, omega_i, w_i))
    index = np.argsort(E0) #Get the index list that would sort E and V
    E0 = E0[index]
    V0 = V0[:,index]
    
    psi0 = V0[:,1]  # Initial state 
    # column 0 -> collective spin pointing along +z
    # column 1 -> collective spin pointing along -z
    
    #######################
    #
    # Parameters of the Krylov
    # basis generation
    #
    #######################

    a = np.zeros(n_krylov) #Lanczos coefficients
    b = np.zeros(n_krylov)
    
    K = np.zeros((dim, n_krylov))  #It stores the Krylov states in the columns
    
    
    #######################
    #
    # Calculating the complexity
    #
    #######################

    
    C = np.zeros((len(omega_f) -1, len(w_f)))  # Stores the time averaged complexity
    
    vect = np.linspace(ti, tf, n_times)
    

    for OO in range(len(omega_f)-1, 0, -1): #run backwards through the values of omega_f
        for ww in range(len(w_f)):
    
            H = H_LMG(chi_f, omega_f[OO], w_f[ww])
            
            ######################
            #
            # Generating the Krylov basis
            #
            ######################
            
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
            
            K, R = np.linalg.qr(K)
            
            
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
                
                
            C[len(omega_f)-1-OO,ww] = (1/tf)*np.trapz(C_t, vect)
            
            # display(Math(r"$\Omega_f = {}$".format(omega_f[OO])))
            # print("(omega_f, w_f) = ({},{})".format(omega_f[OO], w_f[ww]))
            
    return(C, K)

i = 0

fig, ax = plt.subplots()
# ax2 = ax.inset_axes([0.15, 0.35, 0.3, 0.6])
# ax2.tick_params(labelsize=12)

for omega_i in [0.0, 0.1, 0.2, 0.3]:
    j = 100
    omega_f = np.linspace(omega_i, 1, 40)
    w_i = 0
    w_f = [0]
    n_krylov = int(2*j+1)
    ti = 0
    tf = 150
    n_times = 350
    
    # for j in [300]:
        # for omega_i in [0.4]:
            
            # omega_f = np.linspace(omega_i, 2, 30)
            
            # C, K = diagram_complexity(omega_i,
            #                         omega_f,
            #                         w_i,
            #                         w_f,
            #                         n_krylov,
            #                         ti,
            #                         tf,
            #                         n_times)
            
    
    C = np.loadtxt("2DDiagram_Complexity_N{}_h0{}.txt".format(2*j, omega_i))
    
    
    c = ['blue', 'darkorange', 'green', 'red']
    l_style = ['solid', 'dotted', 'dashed', 'dashdot']
    
    # if len(w_f) > 1:
    #     plt.imshow(C/j, cmap='rainbow',
    #                 extent =[w_f.min(), w_f.max(), omega_f.min(), omega_f.max()], 
    #                 aspect='auto')
    #     plt.colorbar(orientation="vertical")
    #     plt.xlabel(r"$\omega$", fontsize=17)
    #     plt.ylabel(r"$\Omega$", fontsize=17)
    #     plt.title(r"$2\overline{C}/N$", loc='right')
    #     plt.tick_params(length = 1)
    #     plt.tight_layout()
    #     # plt.savefig("2D_diagram_TimeAveragedComplexity_N{}.pdf".format(2(j)))
    
    # else:
    C = C[::-1]
    ax.plot(omega_f[1:], C/j, color=c[i], linestyle=l_style[i], label=r"$h_0={}$".format(omega_i))
    ax.axvline((omega_i+1)/2, linestyle=l_style[i], color=c[i], linewidth=1.5)
    
    
    # ax2.plot((omega_f[2:]+omega_f[1:-1])[13:22]/2, np.diff(C)[13:22]/np.diff(omega_f[1:])[13:22], color=c[i], linestyle=l_style[i])
    # ax2.axvline((omega_i+1)/2, linestyle=l_style[i], color=c[i], linewidth=1.5)
    # ax2.set_ylim([0, 4000])
    i += 1


# ax2.axvline((omega_i+1)/2, linestyle='--', color='black')

ax.set_xlabel(r"$h_f$", fontsize=20)
# ax2.set_xlabel(r"$h_f$", fontsize=12)
ax.set_ylabel(r"$\overline{C}/j$", fontsize=20)
# ax2.set_ylabel(r"$d\overline{C}/dh_f$", fontsize=12)
ax.set_ylim([0, 1.05])
plt.text(0.95, 0.05, r'$(a)$', fontsize=22)
ax.legend(frameon=False, fontsize=16)
plt.tight_layout()

plt.savefig("TimeAveragedComplexity_N{}.eps".format(2*j))

        







