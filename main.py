import numpy as np
import sweeper
import math

#to add: multigroups in energy, non diamond differencing 
segments_number= 6000  
angle_number=10     
slab_size= 1 
EnGroup_number= 30      
SI_number=100  #number of internal (source) iterations
PN_order=5 

e_source=2/3 

def external_source(u,x):
    return math.exp(-x)*(4/3-u)

def boundary_condition_1(u):
    return 1
def boundary_condition_2(u):
    return 1/math.exp(1)       #evaluated on negative cosines

b_conditions=[boundary_condition_1,boundary_condition_2]

psi=np.zeros((segments_number,angle_number))
psi_tilde=sweeper.transport_sweep(psi,slab_size,PN_order, b_conditions,external_source)
psi=psi_tilde


for i in range(SI_number):
    psi_tilde=sweeper.transport_sweep(psi_tilde,slab_size,PN_order)
    psi+=psi_tilde
print(psi[0])
