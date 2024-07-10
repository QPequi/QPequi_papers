Como salvar dados em .txt?

    > Importe a biblioteca numpy para trabalhar com arrays numéricos e a biblioteca os para manipulação de arquivos.
    > Defina o intervalo de x=[0,5] com 100 pontos usando a função linspace() do numpy.
    > Calcule o vetor y com a fórmula y = ax^2 + 3x, onde a é um valor qualquer que você pode definir.
    > Salve o vetor y em um arquivo .txt usando a biblioteca os.

Código:

import numpy as np
import os

# Defina o intervalo de x
x = np.linspace(0, 5, 100)

# Defina o valor de a
a = 2

# Calcule o vetor y
y = a*x**2 + 3*x

# Salve o vetor y em um arquivo .txt
np.savetxt(os.path.join(os.getcwd(), 'vetor_y.txt'), y)

    Nesse exemplo, o valor de a foi definido como 2, mas você pode alterá-lo para o valor que desejar. O arquivo .txt será salvo no diretório atual de trabalho com o nome "vetor_y.txt". Note que estou utilizando a função os.path.join para garantir que o arquivo seja salvo no sistema de arquivos correto independente do sistema operacional.

Como usar o arquivo .txt para plotar o gráfico?

    > Importe as bibliotecas numpy e matplotlib.
    > Carregue o vetor y a partir do arquivo .txt usando a função loadtxt() do numpy.
    > Defina o intervalo de x=[0,5] com 100 pontos usando a função linspace() do numpy.
    > Use a função plot() do matplotlib para plotar o gráfico de y em função de x.

Código:

import numpy as np
import matplotlib.pyplot as plt

# Carregue o vetor y a partir do arquivo .txt
y = np.loadtxt('vetor_y.txt')

# Defina o intervalo de x
x = np.linspace(0, 5, 100)

# Plot o gráfico de y em função de x
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Gráfico de y = a*x^2 + 3*x')
plt.show()

    Esse código irá plotar o gráfico de y em função de x, com as respectivas legendas nos eixos x e y e um título para o gráfico. Note que você pode alterar essas legendas e o título para o que desejar.






















