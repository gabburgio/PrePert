import numpy as np
import Xsection_moments
import source_calculator
Sigma_t=1.5  

def transport_sweep(old_psi,slab_size,PN_order, boundary_conditions=None, ext_source=None):

    segments_number=np.shape(old_psi)[0]
    angle_number=np.shape(old_psi)[1]
    points=np.polynomial.legendre.leggauss(angle_number)[0]

    #get cross section moments
    cross_section_moments=Xsection_moments.calculate_Xsection_moments(10,10)

    #get new scattering source
    source=source_calculator.Calculate_source(old_psi,PN_order,cross_section_moments,ext_source)

    #new flux calculation
    new_psi=np.zeros(np.shape(old_psi))
    for i in range(angle_number):        
        if(points[i]>=0):
            if(boundary_conditions is not None):
                psi_minus=boundary_conditions[0](points[i])
            else:
                psi_minus=0
            for j in range(segments_number):
                new_psi[j][i]=(source[j][i]+psi_minus*(2*points[i]*segments_number/slab_size))/(2*points[i]*segments_number/slab_size + Sigma_t)   #calculate new angular flux value
                psi_minus=2*new_psi[j][i]-psi_minus     #get new boundary value (psi_plus) for next spatial segment
        else:
            if(boundary_conditions is not None):            
                psi_plus=boundary_conditions[1](points[i])
            else:
                psi_plus=0
            for j in reversed(range(segments_number)):
                new_psi[j][i]=(source[j][i]-psi_plus*(2*points[i]*segments_number/slab_size))/(-2*points[i]*segments_number/slab_size + Sigma_t)   #the denominator is always positive because the cosine is negative
                psi_plus=2*new_psi[j][i]-psi_plus     #get new boundary value (psi_minus) for next spatial segment
        
    return new_psi



    
    



