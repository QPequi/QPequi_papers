import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
import math as mt
# import os

# This code compute the exact lower bound of entropy production in terms of function s(x)
dim2 = 300
ll = 5
j = ll*100
dim = int(2*j+1)
g = 0.5
vect = np.linspace(0,30,dim2)

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

L = np.zeros(len(vect))  # Loschmidt Echo
LL = np.zeros(len(vect)) # Bures Angle
s_ex_min = np.zeros(dim2)
s_ex = np.zeros((dim2,dim2))
s_lb = np.zeros(dim2)
s_ub = np.zeros(dim2)

#Final hamiltonian

h = 0.8
H = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*h*Jz(j)

e = 10**(-12)

for tt in range(len(vect)):
    U = expm(-1j*H*vect[tt])
    L[tt] = np.abs(np.dot(dag(V0[:,0]), np.dot(U, V0[:,0])))**2
    LL[tt] = (mt.pi/2)*mt.acos(L[tt])#já transformado em x
    #LL[tt] = mt.acos(L[tt])
    s_lb[tt] = (LL[tt]**2)/3
    s_ub[tt] = - np.log(1-LL[tt]*(2/mt.pi)**2)


for xx in range(len(LL)):
      rvec = np.linspace(((2/mt.pi)*LL[xx])+e,(mt.pi/2)-e,dim2)
      for rr in range(len(rvec)):
          s_ex[xx,rr] = (rvec[rr] + LL[xx])*np.log((rvec[rr]+LL[xx])/rvec[rr]) + (1+rvec[rr]-LL[xx])*np.log((1+rvec[rr]-LL[xx])/(1+rvec[rr]))
         
for ii in range(len(LL)):
      s_ex_min[ii] = np.min(s_ex[ii,]) 

plt.plot(vect,s_ex_min,label='exact')
plt.plot(vect,s_lb,label='lower')
plt.plot(vect,s_ub,label='upper')
plt.legend()
plt.savefig('s(x).svg', format='svg')

#np.savetxt(os.path.join(os.getcwd(),'s_exact.txt'), s_exact_min)
