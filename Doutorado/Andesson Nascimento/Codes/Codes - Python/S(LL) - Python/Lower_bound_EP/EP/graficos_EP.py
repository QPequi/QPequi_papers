import numpy as np
import matplotlib.pyplot as plt
import time

# Início da contagem de tempo
start_time = time.time()

# Carregue o vetor y a partir do arquivo .txt

EP_h01 = np.loadtxt('EP_h=0.1_j=500.txt')
EP_h02 = np.loadtxt('EP_h=0.2_j=500.txt')
EP_h04 = np.loadtxt('EP_h=0.4_j=500.txt')
EP_h05 = np.loadtxt('EP_h=0.5_j=500.txt')
EP_h06 = np.loadtxt('EP_h=0.6_j=500.txt')
EP_h08 = np.loadtxt('EP_h=0.8_j=500.txt')
EP_h10 = np.loadtxt('EP_h=1_j=500.txt')

# Defina o intervalo de x
vect = np.linspace(0,30,200)

# Plot o gráfico de y em função de x

plt.plot(vect, EP_h08)
plt.xlabel('t', fontsize=14)
plt.ylabel(r'$\langle \Sigma \rangle$', fontsize=16)
plt.title('h = 0,8')
# plt.legend()
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
plt.savefig('EP_h=0,8_j=500.svg', format='svg')






