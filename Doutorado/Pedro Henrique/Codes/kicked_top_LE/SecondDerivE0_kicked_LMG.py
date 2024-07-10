


import numpy as np
from numpy import linalg as LA
# from scipy.linalg import fractional_matrix_power
# from scipy.linalg import expm
# from scipy.linalg import logm
# from scipy.signal import find_peaks
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 10
dim = np.int(2*j+1)
gamma = 1


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



##################################################
#
# Second derivative of the ground state energy as
# as a function of h
#
##################################################

L = 150  #number of h values

En = np.zeros((dim, L))  #store the ground state energies

h_values = np.linspace(0,1,L)

i = 0
for h in h_values:    
    H = -(1-h)*Jz(j) - h/2/j*np.dot(Jx(j),Jx(j)) 
    
    #Hamiltonian's eigenvalues and eigenvectors
    E, V = LA.eig(H)
    index = np.argsort(E) #Get the index list that would sort E and V
    E = E[index] 
    V = V[:,index]

    En[:,i] = E/j
    i = i+1
    

# # First derivative of E0
# dy = np.diff(E0)
# dx = np.diff(h_values)
# yfirst = dy/dx
# xfirst = 0.5*(h_values[:-1] + h_values[1:])

# # Second derivative of E0
# ddy = np.diff(yfirst)
# ddx = np.diff(xfirst)
# ysecond = ddy/ddx
# xsecond = 0.5*(xfirst[:-1] + xfirst[1:])

for j in range(dim):
    plt.plot(h_values, En[j,:])
    
plt.xlabel(r"$h$", fontsize=17)
plt.ylabel(r"$E_0/j$", fontsize=17)
plt.tight_layout()
plt.savefig("Spectrum_vs_h_LMG_j10.pdf")

# plt.plot(xsecond, (1/j)*ysecond, label=r"$j={}$".format(j))
# plt.legend(fontsize=12, loc=4)
# plt.ylabel(r"$\frac{1}{j}\frac{d^2 E_0}{dh^2}$", fontsize=17)
# plt.xlabel(r"$h$", fontsize=17)
# # plt.xlim([0, 2])
# plt.grid(True)
# plt.tight_layout()
    
# plt.savefig("SecondDerivE0_vs_h.pdf")

