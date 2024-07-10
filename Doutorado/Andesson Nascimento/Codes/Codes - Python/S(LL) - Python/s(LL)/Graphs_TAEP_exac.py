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

DerT_ub = np.gradient(TAEP_ub3,hvec)
DerT_lb = np.gradient(TAEP_lb3,hvec)
DerT_exac = np.gradient(TAEP_exac3,hvec)

# Plotting variation of Time Average Entropy Production
# 
# plt.plot(hvec,DerT_ub)
# plt.plot(hvec,DerT_lb)
# plt.plot(hvec,DerT_exac)
# # plt.ylabel(r'$\frac{d\bar{\langle \Sigma \rangle}}{dh}$', fontsize=14)
# plt.xlabel('h', fontsize=14)
# plt.axvline(0.5, color="red",linestyle='--',label='Dynamical Critical Point')
# plt.legend()

# Ploting Time Average Entropy Production

# plt.plot(hvec,TAEP_ub1, label='j=100')
plt.plot(hvec,TAEP_lb1,label='j=100')
# plt.plot(hvec,TAEP_exac1,label='j=100')
# plt.plot(hvec,TAEP_ub3, label='j=300')
plt.plot(hvec,TAEP_lb3,label='j=300')
# plt.plot(hvec,TAEP_exac3,label='j=300')
# plt.plot(hvec,TAEP_ub5, label='j=500')
plt.plot(hvec,TAEP_lb5,label='j=500')
# plt.plot(hvec,TAEP_exac5,label='j=500')
plt.ylabel(r'$\overline{s\left(\frac{2}{\pi}\mathcal{L} \right)}$', fontsize=14)
plt.xlabel('h', fontsize=14)
plt.axvline(0.5, color="red",linestyle='--')
plt.legend()

# Salva a figura em formato SVG
plt.savefig('TAEP_s(x)_lower.svg', format='svg')

# Salve o vetor y em um arquivo .txt
# np.savetxt(os.path.join(os.getcwd(),'Derivada_j=300.txt'), DerT)














