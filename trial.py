import numpy as np
import scipy
from scipy.special import legendre
weights=np.polynomial.legendre.leggauss(7)[1]

legendre_poly_values=np.zeros((5,7))

#legendre_poly_values[3]=scipy.special.eval_legendre(2,(points))
print(np.sum(weights))