import numpy as np
import scipy
from scipy.special import legendre

#Functional expression of the scattering cross section

def Sigma_s(u):     
    return 1
vec_Sigma_s=np.vectorize(Sigma_s)

#calculate Legendre moments of the scattering cross-section

def calculate_Xsection_moments(quadrature_points_number,number_of_moments):
    cross_section_moments=np.zeros(number_of_moments)
    moment_quadrature_points=np.polynomial.legendre.leggauss(quadrature_points_number)[0]
    moment_quadrature_weights=np.polynomial.legendre.leggauss(quadrature_points_number)[1]

    Sigma_values=vec_Sigma_s(moment_quadrature_points)
    for i in range(number_of_moments):
        Integrand_values=Sigma_values*scipy.special.eval_legendre(i,(moment_quadrature_points))
        cross_section_moments[i]=2*np.pi*np.dot(Integrand_values, moment_quadrature_weights)
    return cross_section_moments
