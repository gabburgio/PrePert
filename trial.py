import numpy as np
import scipy
from scipy.special import legendre
points=np.polynomial.legendre.leggauss(7)[0]

legendre_poly_values=np.zeros((5,7))

legendre_poly_values[3]=scipy.special.eval_legendre(2,(points))
print(legendre_poly_values)