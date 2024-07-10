

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
            # 0 -> collective spin pointing up
            # 1 -> collective spin pointing down
            
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
                
            
            # re-orthogonalization process applying Gram-Schmidt explicitly
            
            K, Q = np.linalg.qr(K)
                
            ######################
            #
            # Calculating IPR
            #
            ######################
            
            vect = np.linspace(ti, tf, n_times)
            IPR_K = np.zeros(len(vect)) #IPR in the Krylov basis
            IPR_E = np.zeros(len(vect)) #IPR in the energy basis
            
            for tt in range(len(vect)):
                
                print(tt)
                U = expm(-1j*H*vect[tt])
                
                for k in range(n_krylov):
                    
                    IPR_E[tt] = IPR_E[tt] + np.abs(expval(V0[:,k], U, psi0))**4
                    IPR_K[tt] = IPR_K[tt] + np.abs(expval(K[:,k], U, psi0))**4
            
            
            return(IPR_E, IPR_K)
 
        

w_i = 0
w_f = 0
ti = 0
tf = 50
n_times = 350
vect = np.linspace(ti, tf, n_times)

i = 0

IPR_E = []
IPR_K = []


for j in [100]:
    for omega_i, omega_f in [(0.1, 0.4), (0.2, 0.4), (0.1, 0.8), (0.2, 0.8)]:
        n_krylov = int(2*j+1)
        
        IPR = IPR_Krylov(j,
                        omega_i,
                        omega_f,
                        w_i,
                        w_f,
                        n_krylov,
                        ti,
                        tf,
                        n_times)
        
        IPR_E.append(IPR[0])
        IPR_K.append(IPR[1])
    
# np.savetxt("Krylov_IPR_LMG_h{}_N{}_d{}.txt".format(omega_f[0], 2*j,d), IPR_t)


#######################
#
# Plots
#
#######################

if omega_f <= 1/2:
    c = ['#003366', '#0080FF', '#99CCFF']  #blue tones
elif omega_f > 1/2:
    c = ['#660000', '#FF0000', '#FF9999']  #red tones



fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

ax1.plot(vect, IPR_E[0], color='#0080FF')
ax1.plot(vect, IPR_K[0], color='#003366', linestyle='--')

#Insets
x1, x2, y1, y2 = 22, 28, 0.01, 0.2  # subregion of the original image
axins = ax1.inset_axes([0.3, 0.55, 0.3, 0.4],
xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])

axins.plot(vect[140:200], IPR_E[0][140:200], color='#0080FF')
axins.plot(vect[140:200], IPR_K[0][140:200], color='#003366', linestyle='--')


ax1.indicate_inset_zoom(axins, edgecolor="black")

ax1.text(44, 0.85, r'$(a)$')
ax1.text(34, 0.62, r'$h_0={}$'.format(0.1), fontsize=15)
ax1.text(34, 0.45, r'$h_f={}$'.format(0.4), fontsize=15)
ax1.set_ylabel(r'$IPR(t)$', fontsize=20)

ax2.plot(vect, IPR_E[1], color='#0080FF')
ax2.plot(vect, IPR_K[1], color='#003366', linestyle='--')

#Inset

axins = ax2.inset_axes([0.3, 0.55, 0.3, 0.4],
xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])

axins.plot(vect[140:200], IPR_E[1][140:200], color='#0080FF')
axins.plot(vect[140:200], IPR_K[1][140:200],  color='#003366', linestyle='--')


ax2.indicate_inset_zoom(axins, edgecolor="black")

ax2.text(45, 0.85, r'$(b)$')
ax2.text(34, 0.62, r'$h_0={}$'.format(0.2), fontsize=15)
ax2.text(34, 0.45, r'$h_f={}$'.format(0.4), fontsize=15)
# ax2.set_ylim([0, 0.25])

ax3.plot(vect, IPR_E[2], color='#FF0000')
ax3.plot(vect, IPR_K[2], color='#660000', linestyle='--')

# ax3.plot(vect, np.ones(len(vect)), color='black', linestyle='--', linewidth=1)
ax3.text(45, 0.85, r'$(c)$')
ax3.text(34, 0.62, r'$h_0={}$'.format(0.1), fontsize=15)
ax3.text(34, 0.45, r'$h_f={}$'.format(0.4), fontsize=15)
ax3.set_xlabel(r'$t$', fontsize=20)
ax3.set_ylabel(r'$IPR(t)$', fontsize=20)

ax4.plot(vect, IPR_E[3], color='#FF0000')
ax4.plot(vect, IPR_K[3], color='#660000', linestyle='--')
# ax4.plot(vect, np.ones(len(vect)), color='black', linestyle='--', linewidth=1)
ax4.text(44, 0.85, r'$(d)$')
ax4.text(34, 0.62, r'$h_0={}$'.format(0.2), fontsize=15)
ax4.text(34, 0.45, r'$h_f={}$'.format(0.8), fontsize=15)
ax4.set_xlabel(r'$t$', fontsize=20)
# ax4.set_ylim([-0.05, 2.01])

plt.tight_layout(pad=0.3)

# for ax in fig.get_axes():
    # ax.label_outer()
     
# plt.savefig("KrylovIPR_vs_t_LMG_N{}.eps".format(2*j))


