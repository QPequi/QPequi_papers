import numpy as np
import matplotlib.pyplot as plt
import time

# Início da contagem de tempo
start_time = time.time()

# Carregue o vetor y a partir do arquivo .txt

# j_10 = np.loadtxt('Derivada_j=10.txt')
j_50 = np.loadtxt('Derivada_j=50.txt')
# j_100 = np.loadtxt('Derivada_j=100.txt')
j_200 = np.loadtxt('Derivada_j=200.txt')
# j_300 = np.loadtxt('Derivada_j=300.txt')
j_500 = np.loadtxt('Derivada_j=500.txt')

# Defina o intervalo de x
hvec = np.linspace(0.1,1.2,120)

# Plot o gráfico de y em função de x
#plt.plot(hvec, j_10,label='j=10')
plt.plot(hvec, j_50,label='j=50')
#plt.plot(hvec, j_100,label='j=100')
plt.plot(hvec, j_200,label='j=200')
# plt.plot(hvec, j_300,label='j=300')
plt.plot(hvec, j_500,label='j=500')
plt.xlabel('h', fontsize=14)
plt.ylabel(r'$\frac{d\bar{\langle \Sigma \rangle}}{dh}$', fontsize=16)
#plt.title('Variação da TAEP')
plt.axvline(0.5, color="red",linestyle='--',label='Dynamical Critical Point')
plt.legend()
# plt.show() #não pode ser usado com o save.svg

# Fim da contagem de tempo
end_time = time.time()

# Tempo total de execução
total_time = end_time - start_time

# Imprime o tempo total em segundos
print(f"Tempo total de execução: {total_time} segundos")

fig = plt.gcf()
fig.set_size_inches(8, 6)
fig.set_dpi(300)


# Salva a figura em formato SVG
plt.savefig('DerTAEP.svg', format='svg')






