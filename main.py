import numpy as np
from CPM_class import Potts
from time import time
import matplotlib.pyplot as plt

def InitLattice(n,m, q):
    
    L = np.zeros((n,m)).astype(int)
    size_cells  = [int(n/(3.0*q)), int(m/(3.0*q))]
    
    for i in range(1,q+1):
        L[i*int(n/3.) + (i-1)*size_cells[0]:i*int(n/3.) + i*size_cells[0],
          i*int(m/3.) + (i-1)*size_cells[1]:i*int(m/3.) + i*size_cells[1]] = int(i)
        
    return L

n,m = 50,50
Q = 2                           # Number of cells (not including medium)
MCS = int(100)
freq_plot = 2

L = InitLattice(n,m, Q)

J = np.ones((2*(Q+1)+2,Q+1))
J[0,0], J[0,1], J[0,2] = 10., 38., 38.
J[1,0], J[1,1], J[1,2] = J[0,1], 20., 11.
J[2,0], J[2,1], J[2,2] = J[0,2], 20., 19.
J[:3,:3] *= n/2

model = Potts(Q, L, J=J)

mat = plt.matshow(model.latt, cmap = plt.get_cmap('plasma', Q+1), vmin = -0.5,
                  vmax = Q+0.5)
plt.axis('off')
plt.draw()
plt.pause(0.01)

start_time = time()

for t in np.arange(MCS):

    model.MonteCarloStep()

    if t%freq_plot==0:
        mat.set_data(model.latt)
        plt.draw()
        plt.pause(0.0001)

print("Execution time: {}s.".format(time()-start_time))

