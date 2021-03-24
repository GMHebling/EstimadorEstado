import matplotlib.pyplot as plt
import scipy.sparse as sp
import numpy as np

i = []
j = []
x = []

#fname = 'matnec_34_alfa0.txt'
#fname = 'matnec_123_alfa0.txt'
fname = 'matnec_342_alfa0.txt'

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
#plt.spy(mat, markersize=0.2)
#plt.show()