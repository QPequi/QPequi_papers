import numpy as np
import matplotlib.pyplot as plt
# import math as mt

d = 1
dimh = 100 # Dimension of h
hvec = np.linspace(0.1,0.8,dimh)

TAEP_ub1 = np.loadtxt('TAEP_ub_j=100.txt')
TAEP_lb1 = np.loadtxt('TAEP_lb_j=100.txt') 
TAEP_exac1 = np.loadtxt('TAEP_exac_j=100.txt')

TAEP_ub3 = np.loadtxt('TAEP_ub_j=300.txt')
TAEP_lb3 = np.loadtxt('TAEP_lb_j=300.txt') 
TAEP_exac3 = np.loadtxt('TAEP_exac_j=300.txt')

TAEP_ub5 = np.loadtxt('TAEP_ub_j=500.txt')
TAEP_lb5 = np.loadtxt('TAEP_lb_j=500.txt') 
TAEP_exac5 = np.loadtxt('TAEP_exac_j=500.txt')


# TAEP variation in terms of quench 

DerT_ub = np.gradient(TAEP_ub5,hvec)
DerT_lb = np.gradient(TAEP_lb5,hvec)
DerT_exac = np.gradient(TAEP_exac5,hvec)

# Plotting variation of Time Average Entropy Production
# 
plt.plot(hvec,DerT_ub,label='Upper Bound, j=500')
# plt.plot(hvec,DerT_lb,label='Lower Bound, j=500')
# plt.plot(hvec,DerT_exac, label='Exact, j=500')
# plt.ylabel(r'$\frac{d\bar{\langle \Sigma \rangle}}{dh}$', fontsize=14)
plt.xlabel('h', fontsize=14)
plt.axvline(0.5, color="red",linestyle='--',label='Dynamical Critical Point')
plt.ylabel(r'$\frac{d}{dh}\left(\overline{s\left(\frac{2}{\pi}\mathcal{L} \right)}\right)$', fontsize=14)
plt.xlabel('h', fontsize=14)
plt.legend()


# Salva a figura em formato SVG
plt.savefig('DerTAEP_s(x)_Upper_j=500.svg', format='svg')

# Salve o vetor y em um arquivo .txt
# np.savetxt(os.path.join(os.getcwd(),'Derivada_j=300.txt'), DerT)














