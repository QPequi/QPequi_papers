import numpy as np
import matplotlib.pyplot as plt

# Carregue o vetor y a partir do arquivo .txt

j_50 = np.loadtxt('TAEP_j=50.txt')
j_100 = np.loadtxt('TAEP_j=100.txt')
j_200 = np.loadtxt('TAEP_j=200.txt')
j_500 = np.loadtxt('TAEP_j=500.txt')

# Defina o intervalo de x
hvec = np.linspace(0.1,1.2,100)

# Plot o gráfico de y em função de x
plt.plot(hvec, j_50,label='j=50')
plt.plot(hvec, j_100,label='j=100')
plt.plot(hvec, j_200,label='j=200')
plt.plot(hvec, j_500,label='j=500')
plt.xlabel('h', fontsize=14)
plt.ylabel(r'$\bar{\langle \Sigma \rangle}$', fontsize=16)

#plt.title('Variação da TAEP')
plt.axvline(0.5, color="red",linestyle='--',label='Dynamical Critical Point')
plt.legend()

fig = plt.gcf()
fig.set_size_inches(8, 6)
fig.set_dpi(300)

# Salva a figura em formato SVG
plt.savefig('TAEP_.svg', format='svg')










