import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
import math as mt
import os
import time

# Início da contagem de tempo
start_time = time.time()

# Calculo da derivada da média temporal de produção de entropia

j = 10
dim = int(2*j+1)
g = 0.5
vect = np.linspace(0,10,100)
hvec = np.linspace(0.1,1.2,120) # Be carreful! h =/= 0 because LL becomes indeterminate

# Defining the Jx, Jy and Jz operators
def J(j, op):
    Jop = np.zeros((dim, dim), dtype=complex) if op == 'y' else np.zeros((dim, dim))
    
    for m1 in np.arange(-j, j+1, 1):
        for m2 in np.arange(-j, j+1, 1):
            p1, p2 = int(m1+j), int(m2+j)
            if m1 == m2+1 and op == 'x':
                Jop[p1,p2] = np.sqrt(j*(j+1)-m2*(m2+1))/2
            elif m1 == m2-1 and op == 'x':
                Jop[p1,p2] = np.sqrt(j*(j+1)-m2*(m2-1))/2
            elif m1 == m2+1 and op == 'y':
                Jop[p1,p2] = 1j*np.sqrt(j*(j+1)-m2*(m2+1))/2
            elif m1 == m2-1 and op == 'y':
                Jop[p1,p2] = -1j*np.sqrt(j*(j+1)-m2*(m2-1))/2
            elif m1 == m2 and op == 'z':
                Jop[p1,p2] = -m2
                
    return Jop
def dag(A):
    return np.transpose(np.conjugate(A))
#####################################################################
# Initial conditions
#####################################################################
#Initial hamiltonian
h0 = 0
H0 = -(1/j)*(np.dot(J(j,'x'),J(j,'x')) + g*np.dot(J(j,'y'),J(j,'y'))) - 2*h0*J(j,'z') 
#Initial density matrix of the system
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index]
V0 = V0[:,index]

#####################################################################
# TAEP
#####################################################################

L = np.zeros(len(vect))  # Loschmidt Echo
LL = np.zeros(len(vect)) # Bures Angle
EP = np.zeros(len(vect)) # Minimal Entropy Production - 1 term
TAEP = np.zeros(len(hvec)) #  Time Average - Entropy production
dTAEP = np.zeros(len(vect))

for h in range(len(hvec)):
#Final hamiltonian
    H = -(1/j)*(np.dot(J(j,'x'),J(j,'x')) + g*np.dot(J(j,'y'),J(j,'y'))) - 2*hvec[h]*J(j,'z')
    for tt in range(len(vect)):
        U = expm(-1j*H*vect[tt])
        L[tt] = round(np.abs(np.dot(dag(V0[:,0]), np.dot(U, V0[:,0])))**2,10) 
        if L[tt] <= 0:
            L[tt] = 1e-16
        LL[tt] = mt.acos(L[tt])
        EP[tt] = 8/(mt.pi**2)*LL[tt]**2
    TAEP[h] = (1/np.max(vect))*np.trapz(EP, vect)
    
# TAEP variation in terms of quench 

DerT = np.gradient(TAEP,hvec)

# Plotting Time Average Entropy Production

plt.plot(hvec,DerT)
plt.ylabel(r'$\frac{d\bar{\langle \Sigma \rangle}}{dh}$', fontsize=14)
plt.xlabel('h', fontsize=14)
plt.axvline(0.5, color="red",linestyle='--',label='Dynamical Critical Point')
plt.legend()

# obter o caminho do diretório local do Dropbox
db_path = "/home/anderson/Dropbox/QPequi/Doutorado/Andesson Nascimento/Dados_Python"

# Salve o vetor y em um arquivo .txt
np.savetxt(os.path.join(os.getcwd(),'Derivada_j=10.txt'), DerT)



# for j in range(1, 11):
#     # Seu código aqui
#     nome_arquivo = "j={}.txt".format(j)
#     np.savetxt(os.path.join(os.getcwd(), nome_arquivo), DerT)

# Fim da contagem de tempo
end_time = time.time()

# Tempo total de execução
total_time = end_time - start_time

# Imprime o tempo total em segundos
print(f"Tempo total de execução: {total_time} segundos")




