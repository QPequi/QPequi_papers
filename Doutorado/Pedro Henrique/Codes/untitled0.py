



import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import expm
# from scipy import integrate
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


# dim = int(2*j+1)

omega_i = 0
omega_f = np.linspace(0,1,80)
w_i = 0
w_f = [0]
# n_krylov = dim
ti = 0
tf = 250
n_times = 450
    
    
M25 = np.loadtxt("TimeAveragedMagnetization_w0_N50.txt")
M50 = np.loadtxt("TimeAveragedMagnetization_w0_N100.txt")
M100 = np.loadtxt("TimeAveragedMagnetization_w0_N200.txt")

vect = np.linspace(ti,tf,n_times)

c = ['#000000', '#606060', '#A0A0A0']

i=1
for j in [25,50,100]:
    M = np.loadtxt("TimeAveragedMagnetization_w0_N{}.txt".format(2*j))
    plt.plot(omega_f, M[::-1]/j, color=c[i-1], label=r"$N={}$".format(2*j))
    plt.xlabel(r"$h$", fontsize=20)
    if i==1:
        plt.ylabel(r"$\overline{S_z}/j$", fontsize=20)
    plt.legend(frameon=False, fontsize=16)
               # loc=(0.55, 0))
    i += 1


y = -(np.linspace(0,0.5,30) -np.ones(30))**9 -np.ones(30)

plt.plot(np.linspace(0.5, 1, 5), np.zeros(5), linestyle=(5,(10,3)), color='black',linewidth=0.75)
# plt.plot(np.linspace(0, 0.5, 30), np.tanh(9*(np.linspace(0,0.5,30)-0.5*np.ones(30))),linestyle=(5,(10,3)), color='black', linewidth=0.75)
plt.plot(np.linspace(0, 0.5, 30), y[::-1], linestyle=(5,(10,3)), color='black', linewidth=0.75)
# plt.arrow(x=0.35, y=-0.9, dx=0.1, dy=0, width=.02, facecolor='grey', edgecolor='none')
plt.tight_layout()

plt.savefig("TimeAveragedMagnetization_vs_h.pdf")











