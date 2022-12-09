import numpy as np
import sweeper

#to add: multigroups in energy, non uniform scattering, non diamond differencing, take source average on segments, Q as a function
a= 10   #slab size in cm
J= 30  #number of spatial segments
N= 5    #number of discrete ordinates
G= 30   #number of energy groups
Q=1     #external source
M=200  #number of internal (source) iterations

def boundary_condition_1(u):
    return u
def boundary_condition_2(u):
    return -u
boundary_conditions=[boundary_condition_1,boundary_condition_2]
def Source(x):
    return 1

psi=np.zeros((J,N))
psi_tilde=sweeper.transport_sweep(psi,J,a,True,boundary_conditions,Q)
psi=psi_tilde

for i in range(M):
    psi_tilde=sweeper.transport_sweep(psi_tilde,J,a,False)
    psi+=psi_tilde
print(psi)
print("nne")