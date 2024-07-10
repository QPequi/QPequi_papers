import numpy as np
import matplotlib.pyplot as plt
import math as mt

up = np.loadtxt('s_ub.txt')
lower = np.loadtxt('s_lb.txt') 
exact = np.loadtxt('s_exact.txt')
vec = np.linspace(0, mt.pi/2,10)

plt.plot(vec, up,label='upper bound')
plt.plot(vec, lower,label='lower bound')
plt.plot(vec, exact,label='s(x)')
plt.xlabel(r'$\mathcal{L}$', fontsize=14)
plt.ylabel(r'$s(\frac{2}{\pi}\mathcal{L})$', fontsize=16)
#plt.title('Variação da TAEP')
plt.legend()




# plt.savefig('DerTAEP.svg', format='svg')
