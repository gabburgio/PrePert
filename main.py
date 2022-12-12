import numpy as np
import sweeper

#to add: multigroups in energy, non diamond differencing, Q as a function
slab_size= 10   
segments_number= 30  
angle_number= 5    
EnGroup_number= 30   
e_source=1     
SI_number=50  #number of internal (source) iterations
PN_order=5

def boundary_condition_1(u):
    return 1
def boundary_condition_2(u):
    return 1       #evaluated on negative cosines

b_conditions=[boundary_condition_1,boundary_condition_2]

psi=np.zeros((segments_number,angle_number))
psi_tilde=sweeper.transport_sweep(psi,slab_size,PN_order, b_conditions,e_source)
psi=psi_tilde


for i in range(SI_number):
    psi_tilde=sweeper.transport_sweep(psi_tilde,slab_size,PN_order)
    psi+=psi_tilde
print(psi)
