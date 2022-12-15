import numpy as np
import scipy
from scipy.special import legendre
weights=np.polynomial.legendre.leggauss(10)[1]

legendre_poly_values=np.zeros((5,10))

a=np.full((1, 5), 7)
print(a)
#print(sum(weights))