import numpy as np
#from main import J,a

Sigma_s=1   #constant scattering cross-sections
Sigma_t=1.5  

def transport_sweep(boundary_conditions, old_psi, ext_source,J,a):

    #source calculation
    source=np.zeros(np.shape(old_psi)[0])     #number of rows (space segments) in the source
    weights=np.polynomial.legendre.leggauss(np.shape(old_psi)[1])[1]
    for i in range(np.size(source)):
        source[i]=0.5*(Sigma_s*np.dot(weights,old_psi[i])+ext_source)     #source = 1/2(Sigma_s*flux+Q), weighted sum is a scalar product

    #new flux calculation
    points=np.polynomial.legendre.leggauss(np.shape(old_psi)[1])[0]
    new_psi=np.zeros(np.shape(old_psi))
    for i in range(np.shape(old_psi)[1]):        #number of columns (discrete ordinates)
        if(points[i]>=0):
            psi_minus=boundary_conditions[0](points[i])
            for j in range(np.shape(old_psi)[0]):
                new_psi[j][i]=(source[j]+psi_minus*(2*points[i]*J/a))/(2*points[i]*J/a + Sigma_t)   #calculate new angular flux value
                psi_minus=2*new_psi[j][i]-psi_minus     #get new boundary value (psi_plus) for next spatial segment
        else:
            psi_plus=boundary_conditions[1](points[i])
            for j in reversed(range(np.shape(old_psi)[0])):
                new_psi[j][i]=(source[j]-psi_plus*(2*points[i]*J/a))/(-2*points[i]*J/a + Sigma_t)   #calculate new angular flux value
                psi_plus=2*new_psi[j][i]-psi_plus     #get new boundary value (psi_minus) for next spatial segment
        
    return new_psi



    
    



