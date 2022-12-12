import numpy as np
import scipy
from scipy.special import legendre
import Xsection_moments

#get cross section moments
cross_section_moments=Xsection_moments.calculate_Xsection_moments(10,10)

def Calculate_source(old_psi, PN_order, ext_source=None):

    segments_number=np.shape(old_psi)[0]
    angle_number=np.shape(old_psi)[1]
    weights=np.polynomial.legendre.leggauss(angle_number)[1]
    points=np.polynomial.legendre.leggauss(angle_number)[0]

    #source calculation
    source=np.zeros((segments_number,angle_number))

    #old flux moments calculation
    legendre_poly_values=np.zeros((PN_order,angle_number))
    for i in range(PN_order):
        legendre_poly_values[i]=scipy.special.eval_legendre(i,(points))
    old_flux_moments=np.zeros((segments_number, PN_order))
    for i in range(segments_number):
        for j in range(PN_order):
            old_flux_moments[i][j]=np.dot(old_psi[i]*legendre_poly_values[j],weights)
    
    
    #scattering source calculation
    for i in range(segments_number):
        for j in range(angle_number):
            for k in range(PN_order):
                source[i][j]+=(k+0.5)*cross_section_moments[k]*old_flux_moments[i][k]*legendre_poly_values[k][j]
    

    #adding external source
    if(ext_source is not None):
        for i in range(segments_number):
            source[i]+=ext_source*weights
            

    return source
