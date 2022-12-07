import numpy as np
import sweeper

#to add: multigroups in energy, non uniform scattering, non diamond differencing, take source average on segments, Q as a function
a= 10   #slab size in cm
J= 500  #number of spatial segments
N= 5    #number of discrete ordinates
G= 30   #number of energy groups
Q=1     #external source
M=200  #number of internal (source) iterations

def boundary_condition_1(u):
    return 0
def boundary_condition_2(u):
    return 0
boundary_conditions=[boundary_condition_1,boundary_condition_2]
def Source(x):
    return 1

psi=np.zeros((J,N))
psi_tilde=sweeper.transport_sweep(boundary_conditions, psi,Q,N,a)
psi=psi_tilde

for i in range(M):
    psi_tilde=sweeper.transport_sweep(boundary_conditions, psi_tilde,0,N,a)
    psi+=psi_tilde
print(psi)
