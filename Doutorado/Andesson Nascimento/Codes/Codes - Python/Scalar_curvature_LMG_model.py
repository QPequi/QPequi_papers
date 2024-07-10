import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

# sy.init_printing(use_unicode=False, wrap_line=False)
# plt.rcParams["font.size"] = 15
# plt.rcParams["text.usetex"] = True

# Definição dos parâmetros
j = 120
dim = int(2*j+1)
O = 1

T = np.linspace(0,6,10)
EY = np.linspace(0,6,10)

#Defining the Jx, Jy and Jz in z basis
def Jx(j):
    Jx = np.zeros((dim, dim))
    
    aut1 = 0
    aut2 = 0
    
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

# Definição das derivadas
dH_Ox =  Jx(j)
dH_Ey = (1/j)*np.dot(Jy(j),Jy(j))

#Initial hamiltonian
Q11 = np.zeros((len(T),len(EY)))
Q12 = np.zeros((len(T),len(EY)))
Q22 = np.zeros((len(T),len(EY)))
raz_Q = np.zeros((len(T),len(EY)))
R = np.zeros((len(T),len(EY)))
inv_sq = np.zeros((len(T),len(EY)))
A = np.zeros((len(T),len(EY)))
B = np.zeros((len(T),len(EY)))

# Derivações
d1_inv = np.zeros((len(T),len(EY)))
d1_raz = np.zeros((len(T),len(EY)))
d1_g11 = np.zeros((len(T),len(EY)))
d1_g22 = np.zeros((len(T),len(EY)))
d1_g12 = np.zeros((len(T),len(EY)))
d2_inv = np.zeros((len(T),len(EY)))
d2_raz = np.zeros((len(T),len(EY)))
d2_g11 = np.zeros((len(T),len(EY)))
D12_g11 = np.zeros((len(T),len(EY)))
D11_g22 = np.zeros((len(T),len(EY)))
D22_g11 = np.zeros((len(T),len(EY)))
D21_g11 = np.zeros((len(T),len(EY)))
D21_g12 = np.zeros((len(T),len(EY)))

k = 0
for Ey in EY:
    i = 0
    for Ox in T:
        H = O*Jz(j) + Ox*Jx(j)+ (Ey/j)*np.dot(Jy(j),Jy(j))
        #Initial hamiltonian's eigenvalues and eigenvectors
        E, V = LA.eig(H)
        index = np.argsort(E) #Get the index list that would sort E and V
        E = E[index] 
        V = V[:,index]
        for m in range(dim):  # n=0 (see article)
            if m != 0:
                    # Elementos da Métrica
                Q11[i,k] = Q11[i,k] + \
                (np.dot(V[:,0],np.dot(dH_Ox, V[:,m]))*np.dot(V[:,m],np.dot(dH_Ox,V[:,0])))/(E[m]-E[0])**2
                Q12[i,k] = Q12[i,k] + \
                (np.dot(V[:,0], np.dot(dH_Ox, V[:,m]))*np.dot(V[:,m],np.dot(dH_Ey,V[:,0])))/(E[m]-E[0])**2
                Q22[i,k] = Q22[i,k] + \
                (np.dot(V[:,0], np.dot(dH_Ey, V[:,m]))*np.dot(V[:,m],np.dot(dH_Ey,V[:,0])))/(E[m]-E[0])**2
                    # Inverso da raiz do determinante da Métrica
                inv_sq[i,k] = 1/np.sqrt(Q11[i,k] + Q22[i,k] - 2*Q12[i,k])
                    # Razão Q12/Q11
                raz_Q[i,k] = Q12[i,k]/Q11[i,k]
        i = i+1
    k = k+1

kk=0
for kk in range(len(EY)):
    ii = 0
    for ii in range(len(T)):
        # Termos de R
        # Derivadas
        d1_inv[kk] = np.gradient(inv_sq[:,kk])
        d1_raz[kk] = np.gradient(raz_Q[:,kk])
        d1_g11[kk] = np.gradient(Q11[:,kk])
        d1_g22[kk] = np.gradient(Q22[:,kk])
        d1_g12[kk] = np.gradient(Q12[:,kk])
        d2_inv[ii] = np.gradient(inv_sq[ii,:])
        d2_raz[ii] = np.gradient(raz_Q[ii,:])
        d2_g11[ii] = np.gradient(Q11[ii,:])
        D12_g11[kk] = np.gradient(d2_g11[:,kk])
        D11_g22[kk] = np.gradient(d1_g22[:,kk])
        D21_g12[ii] = np.gradient(d1_g12[ii,:])
        D22_g11[ii] = np.gradient(d2_g11[ii,:])
        D21_g11[ii] = np.gradient(d1_g11[ii,:])
        
        # Elementos
        A[ii,kk] = raz_Q[ii,kk]*d1_inv[ii,kk]*d2_g11[ii,kk] - \
            d1_inv[ii,kk]*d1_g22[ii,kk] + \
            inv_sq[ii,kk]*d1_raz[ii,kk]*d2_g11[ii,kk] + \
            inv_sq[ii,kk]*raz_Q[ii,kk]*D12_g11[ii,kk] - \
            inv_sq[ii,kk]*D11_g22[ii,kk]  
        B[ii,kk] = 2*d2_inv[ii,kk]*d1_g12[ii,kk] - \
            d2_inv[ii,kk]*d2_g11[ii,kk] - \
            raz_Q[ii,kk]*d2_inv[ii,kk]*d1_g11[ii,kk] + \
            2*inv_sq[ii,kk]*D21_g12[ii,kk] - \
            inv_sq[ii,kk]*D22_g11[ii,kk] - \
            inv_sq[ii,kk]*d2_raz[ii,kk]*d1_g11[ii,kk] - \
            inv_sq[ii,kk]*raz_Q[ii,kk]*D21_g11[ii,kk]                
        R[ii,kk] = inv_sq[ii,kk]*(A[ii,kk]+B[ii,kk])
        ii = ii+1
    kk = kk+1

plt.plot(T,R)

# #plot 3D

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.contour3D(T, EY, Q11, 150)
# ax.contour3D(T, EY, Q12, 150)
# ax.contour3D(T, EY, R, 150)


# fig, ax = plt.subplots(2,2,figsize =(15,10))
# ax[0][0].plot(T,Q11,label=Q11) # 1º grafico
# ax[0][1].plot(T,Q12,label=Q12) # 2º grafico
# ax[1][0].plot(T,Q22,label=Q22) # 3º grafico
# ax[1][1].plot(T,R) # 4º grafico
# plt.legend()

