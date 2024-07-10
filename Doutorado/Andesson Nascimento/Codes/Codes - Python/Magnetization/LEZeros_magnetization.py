
##########################################
#
# This code relates the LE zeros with the
# average magnetization zeros
#
###########################################3



import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
from scipy.signal import find_peaks


plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = True


j = 20
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


def dag(A):
    return np.transpose(np.conjugate(A))


def find_zeros(A, T):
    Dif = np.diff(np.sign(A))
    zerotime = []
    for i in range(len(Dif)):
        if A[i] != 0:
            zerotime.append((T[i+1]+T[i])/2)

    return zerotime



J = 1
g = 0


# Initial hamiltonian (critical point h = 1)
h0 = 0
H0 = -J/j*(np.dot(Jx(j), Jx(j)) + g*np.dot(Jy(j), Jy(j))) - h0*Jz(j) + (1+g)*J/2*np.identity(dim)


#Initial hamiltonian's eigenvalues and eigenvectors
E0, V0 = LA.eig(H0)
index = np.argsort(E0) #Get the index list that would sort E and V
E0 = E0[index]
V0 = V0[:,index]


#Quenched hamiltonian (critical point h = 1)
h = 0.6
H = -J/j*(np.dot(Jx(j), Jx(j)) + g*np.dot(Jy(j), Jy(j))) - h*Jz(j) + (1+g)*J/2*np.identity(dim)


#########################################
#
# LE zeros and magnetization
#
#########################################

vect = np.linspace(0, 50, 150)

M = np.zeros(len(vect)) #magnetization vector
R= np.zeros(len(vect)) #Loschmidt echo vector


for tt in range(len(vect)):
    U = expm(-1j*vect[tt]*H)

    psit = np.dot(U, V0[:,0])

    M[tt] = np.dot(np.conjugate(psit), np.dot(Jx(j), psit)) #magnetization at time vect[tt]

    R[tt] = -(1/2/j)*np.log(np.abs(np.dot(np.conjugate(V0[:,0]), psit))**2) #Rate function at time vect[tt]



t_c, _ = find_peaks(R, height=0) #find the peaks of the rate function
t_m = find_zeros(M, vect) #find the zeros of the magnetization
L = min(len(t_c), len(t_m))


plt.plot(vect[t_c[:L]], t_m[:L], label=r"$N = {}, h_f={}$".format(2*j, h))
# plt.plot(vect, M, label=r"$j = {}, h_f={}$".format(j, h))
# plt.plot(vect[t_c], R[t_c], 'x', c="red")
plt.legend(fontsize = 12)
plt.xlabel(r"$t_M$", fontsize=18)
plt.ylabel(r"$t_c$", fontsize=18)
plt.tight_layout()
plt.show()

plt.savefig("CriticalTimes_vs_ZeroMagnetization_.pdf")


# plt.plot(vect, M/max(M), label=r"$N={}, h_f={}$".format(2*j, h))
# plt.legend(fontsize=12)
# plt.xlabel(r"$t$", fontsize=18)
# plt.ylabel(r"$\langle S_x \rangle (t)$", fontsize=18)
# plt.tight_layout()

# plt.savefig("Magnetization_h0.6_N200.pdf")
