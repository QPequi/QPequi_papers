import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

# sy.init_printing(use_unicode=False, wrap_line=False)
# plt.rcParams["font.size"] = 15
# plt.rcParams["text.usetex"] = True

# Definição dos parâmetros
j = 120
dim = int(2*j+1)
Ey = 0.01
O = 1

T = np.linspace(0,6,100)
#E = np.linspace(0,5,100)

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

Q11 = np.zeros(len(T))
Q12 = np.zeros(len(T))
Q22 = np.zeros(len(T))
sq_detg = np.zeros(len(T))

i=0
for Ox in T:
    H = O*Jz(j) + Ox*Jx(j)+ (Ey/j)*np.dot(Jy(j),Jy(j))
    #Initial hamiltonian's eigenvalues and eigenvectors
    E, V = LA.eig(H)
    index = np.argsort(E) #Get the index list that would sort E and V
    E = E[index] 
    V = V[:,index]
    for m in range(dim):  # n=0 (see article)
        if m != 0:
            Q11[i] = Q11[i] + (np.dot(V[:,0], np.dot(dH_Ox, V[:,m]))*np.dot(V[:,m],np.dot(dH_Ox,V[:,0])))/(E[m]-E[0])**2
            Q12[i] = Q12[i] + (np.dot(V[:,0], np.dot(dH_Ox, V[:,m]))*np.dot(V[:,m],np.dot(dH_Ey,V[:,0])))/(E[m]-E[0])**2
            Q22[i] = Q22[i] + (np.dot(V[:,0], np.dot(dH_Ey, V[:,m]))*np.dot(V[:,m],np.dot(dH_Ey,V[:,0])))/(E[m]-E[0])**2
            sq_detg[i] = Q11[i] + Q22[i] - 2*Q12[i]
    i = i+1

## Resultados analíticos

g_11 = np.zeros(len(T))
g_12 = np.zeros(len(T))
g_22 = np.zeros(len(T))

k=0
for Ox in T:
    g_11[k] = j/(2*(Ox**2+1)**(7/4)*np.sqrt((np.sqrt(Ox**2+1)+2*Ey))) + (Ey**2*Ox**2)/(8*(Ox**2+1)**2*(np.sqrt(Ox**2+1)+2*Ey)**2)
    g_12[k] = -Ey*Ox/(8*(Ox**2+1)*(np.sqrt(Ox**2+1)+2*Ey)**2)
    g_22[k] = 1/(8*(np.sqrt(Ox**2+1)+2*Ey)**2)
    k = k+1

R_an = -4

### Definições para escrever R
# diff_q11 = np.gradient(Q11)
# diff_q22 = np.gradient(Q22)
# diff_q12 = np.gradient(Q12)
# R = np.zeros(len(T))
# for l in range(len(T)):
#     R[l] = 1/sq_detg[l]


## Gráficos
fig, ax = plt.subplots(2,2,figsize =(15,10))
ax[0][0].plot(T,Q11,marker='x',label='Numérico') # 1º gráfico
ax[0][0].plot(T,g_11,color='r',label='Analítico') 
ax[0][0].legend()
ax[0][1].plot(T,Q12,marker='x',label='Numérico') # 2º gráfico
ax[0][1].plot(T,g_12,color='r',label='Analítico')
ax[0][1].legend()
ax[1][0].plot(T,Q22,marker='x',label='Numérico') # 3º gráfico
ax[1][0].plot(T,g_22,color='r',label='Analítico')
ax[1][0].legend()
# ax[1][1].plot(T,R_an,marker='x',label='Numérico') # 4º grafico
ax[1][1].plot(T,R_an,color='r',label='Analítico') 
ax[1][1].legend()
plt.legend()













