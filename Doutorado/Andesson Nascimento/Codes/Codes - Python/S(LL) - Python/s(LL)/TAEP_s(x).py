import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
import math as mt
import os

dim = 50 # Dimension of time vector
dimh = 100 # Dimension of h

vect = np.linspace(0,12,dim)
hvec = np.linspace(0.1,0.8,dimh)
TAEP_lb = np.zeros(len(hvec)) #  Time Average - lower bound
TAEP_ub = np.zeros(len(hvec)) #  Time Average - upper bound
TAEP_exac = np.zeros(len(hvec)) #  Time Average - exac

s_lb = np.zeros(dim)
s_ub = np.zeros(dim)
s_exac = np.zeros((dim,dim))
s_exac_min = np.zeros(dim)

# Computing Bures Angle

# This parth of the code computes the Bures Angle (LL) to initial hamiltonian  H0 = -(1/j)*(np.dot(Jx(j),Jx(j))
# + g*np.dot(Jy(j),Jy(j))) - 2*h0*Jz(j) with J = 1 (and dynamical critical point h = 1)

ll = 3
j = ll*100
dimm = int(2*j+1)
g = 0.5

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

for h in range(len(hvec)):
    #Final hamiltonian
    H = -(1/j)*(np.dot(Jx(j),Jx(j)) + g*np.dot(Jy(j),Jy(j))) - 2*hvec[h]*Jz(j)
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
        if s_ub[tt] > 10^2:
            s_ub[tt] = 30 #Escolhi esse valor pois a função s(x) assume valores próximos de 30 antes de se tornar infinito
    e = 1*10**-12 
            
    for tt in range(len(vect)):
        rvec = np.linspace(((2/mt.pi)*LL[tt])+e,1-e,dim)
        for rr in range(len(rvec)):
            s_exac[tt,rr] = ((rvec[rr]-(2/mt.pi)*LL[tt])*np.log((rvec[rr]-(2/mt.pi)*LL[tt])/rvec[rr]) +\
                              + (1-rvec[rr]+(2/mt.pi)*LL[tt])*np.log((1-rvec[rr]+(2/mt.pi)*LL[tt])/(1-rvec[rr])))                   
            
                
    for ii in range(len(LL)):
        s_exac_min[ii] = np.min(s_exac[ii,])


    # Valor que você deseja substituir os NaNs
    valor_substituto = 30

    # Loop for para percorrer a array
    for ii in range(len(s_exac_min)):
        # Verifica se o valor é NaN
        if np.isnan(s_exac_min[ii]):
            # Substitui o NaN pelo valor definido
            s_exac_min[ii] = valor_substituto


    TAEP_exac[h] = (1/np.max(vect))*np.trapz(s_exac_min, vect) # média de s_exac
    TAEP_lb[h] = (1/np.max(vect))*np.trapz(s_lb, vect) # média de s_lb
    TAEP_ub[h] = (1/np.max(vect))*np.trapz(s_ub, vect) # média de s_ub


plt.plot(hvec,TAEP_ub, label='upper')
plt.plot(hvec,TAEP_lb,label='lower',marker='o')
plt.plot(hvec,TAEP_exac,label='s(x)',marker='x')
plt.axvline(0.5, color="red",linestyle='--')
plt.legend()

np.savetxt(os.path.join(os.getcwd(),'TAEP_exac_j={}.txt'.format(j)), TAEP_exac)
np.savetxt(os.path.join(os.getcwd(),'TAEP_lb_j={}.txt'.format(j)), TAEP_lb)
np.savetxt(os.path.join(os.getcwd(),'TAEP_ub_j={}.txt'.format(j)), TAEP_ub)

# Notas sobre o código:
#   Neste código estou calculando a média temporal para a função s(x), tanto para os limites superior e inferior 
#   quanto para o valor exato, já usando x = (2/pi)*LL, com LL sendo o Ângulo de Bures.
#   O limite inferior é calculado sem problemas. Não há divergências no cálculo de s(x), não gerando divergência para
#   média temporal.

#   O limite superior e o valor exato de s(x) apresentam divergências, o que gera divergências no cálculo da média temporal

#   Para o limite superior, a média temporal para valores elevados de h torna-se infinita devido aos infinitos que surgem na função s(x)
#   para valores os valores de h correspondentes.
#   Para o valor exato, a média temporal de s(x) apresenta valores nan sempre que s(x) apresentar valores nan.

#   Análise dos valores divergentes:
    # Limite superior: Sempre que LL = 1.57079633 (pi/2), o valor de s(x) é infinito (temos log(0)). 
    # Nas proximidades de infinitos, o valor máximo é em torno de 30

    # Valor exato: Sempre que LL = 1.57079633 (pi/2), o valor de s(x) é nan
    # Nas proximidades de nan, os valores de s(x) ficam em torno de 30 (para j=100)

































