import numpy as np
import matplotlib.pyplot as plt

j_100 = np.loadtxt('Derivada_j=100.txt')
j_200 = np.loadtxt('Derivada_j=200.txt')
j_500 = np.loadtxt('Derivada_j=500.txt')

j_100_part = j_100[33:53]
j_200_part = j_200[33:53]
j_500_part = j_500[33:53]


# Defina o intervalo de x
hvec = np.linspace(0.1,1.2,120)

hvec_p = hvec[33:53]

# Plot o gráfico de y em função de x
plt.plot(hvec_p, j_100_part,label='j=100')
plt.plot(hvec_p, j_200_part,label='j=200')
plt.plot(hvec_p, j_500_part,label='j=500')
plt.xlabel('h', fontsize=14)
plt.ylabel(r'$\frac{d\bar{\langle \Sigma \rangle}}{dh}$', fontsize=16)
plt.axvline(0.5, color="red",linestyle='--',label='Dynamical Critical Point')
plt.legend()

fig = plt.gcf()
fig.set_size_inches(8, 6)
fig.set_dpi(300)


# Salva a figura em formato SVG
plt.savefig('DerTAEP_[0.4-0.6].svg', format='svg')






