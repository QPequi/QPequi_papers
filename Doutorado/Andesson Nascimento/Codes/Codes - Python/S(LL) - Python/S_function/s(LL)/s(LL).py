import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
import math as mt

# This code compute the function s(x) (lower, upper and exact) according paper J. Math. Phys. 46, 102104 (2005).

dim = 50 # Dimension of time vector

s_lb = np.zeros(dim)
s_ub = np.zeros(dim)
s_exac = np.zeros((dim,dim))
s_exac_min = np.zeros(dim)

# Computing Bures Angle

# This parth of the code computes the Bures Angle (LL) to initial hamiltonian  H0 = -(1/j)*(np.dot(Jx(j),Jx(j))
# + g*np.dot(Jy(j),Jy(j))) - 2*h0*Jz(j) with J = 1 (and dynamical critical point h = 1)

ll = 2
j = ll*100
dimm = int(2*j+1)
g = 0.5
vect = np.linspace(0,20,dim)

#Defining the Jx, Jy and Jz in z basis
def Jx(j):
    Jx = np.zeros((dimm, dimm))
    
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
    Jy = np.zeros((dimm, dimm),dtype=complex)
    
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
    Jz = np.zeros((dimm, dimm))
    
    for m in np.arange(-j, j+1, 1):
        p = int(m+j)
        Jz[p,p] = -m
    return Jz
def dag(A):
    return np.transpose(np.conjugate(A))

#Initial hamiltonian

h0 = 0
H0 = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h0*Jz(j)

#Initial density matrix of the system
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index]
V0 = V0[:,index]

#####################################################################
# Bures Angle - LL
#####################################################################

L = np.zeros(len(vect))  # Loschmidt Echo
LL = np.zeros(len(vect)) # Bures Angle

# for h in range(len(hvec)):
#Final hamiltonian
h = 0.5

H = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h*Jz(j)
for tt in range(len(vect)):
    U = expm(-1j*H*vect[tt])
    L[tt] = np.abs(np.dot(dag(V0[:,0]), np.dot(U, V0[:,0])))**2
    LL[tt] = mt.acos(L[tt])
    if mt.acos(L[tt]) > (mt.pi/2):
        LL[tt] = (mt.pi/2)

# This parth of the code computing lower bound, upper bound e exact value of s(x).

for tt in range(len(vect)):
    s_lb[tt] = (8/mt.pi**2)*LL[tt]**2
    s_ub[tt] = - np.log(1-(2/mt.pi)*LL[tt])
    
e = 1*10**-12 
            
for tt in range(len(vect)):
    rvec = np.linspace(((2/mt.pi)*LL[tt])+e,1-e,dim)
    for rr in range(len(rvec)):
        s_exac[tt,rr] = ((rvec[rr]-(2/mt.pi)*LL[tt])*np.log((rvec[rr]-(2/mt.pi)*LL[tt])/rvec[rr]) +\
                        + (1-rvec[rr]+(2/mt.pi)*LL[tt])*np.log((1-rvec[rr]+(2/mt.pi)*LL[tt])/(1-rvec[rr])))
            
                    
for ii in range(len(LL)):
    s_exac_min[ii] = np.min(s_exac[ii,])
    # if s_exac_min[ii] == mt.nan:
    #     s_exac_min[ii] = 1    
                
# # plt.plot(vect,s_lb,label='lower')
# plt.plot(vect,s_ub,label='upper',marker='x')
# plt.plot(vect,s_exac_min,label='exact')
# # plt.ylim(0, 20)
# # plt.xlim(0, 10)
# plt.legend()



























