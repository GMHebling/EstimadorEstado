import matplotlib.pyplot as plt
import scipy.sparse as sp
import numpy as np

i = []
j = []
x = []

#fname = 'matrizH_34.txt'
#fname = 'matrizH_123.txt'
#fname = 'matrizH_342.txt'

with open(fname) as f:
    for line in f:
        aux = line.split(',')
        i.append(int(aux[0]))
        j.append(int(aux[1]))
        x.append(float(aux[2]))

i = np.array(i)
j = np.array(j)
x = np.array(x)

mat = sp.coo_matrix((x, (i,j)))
cond = np.linalg.cond(mat.todense())
print(cond)
plt.spy(mat, markersize=0.2)
plt.show()