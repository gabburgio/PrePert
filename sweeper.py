import numpy as np
#from main import J,a

Sigma_s=1   #constant scattering cross-sections
Sigma_t=1.5  

def transport_sweep(old_psi,J,a,zeroth, boundary_conditions=None, ext_source=None):

    #source calculation
    source=np.zeros(np.shape(old_psi)[0])     #number of rows (space segments) in the source
    weights=np.polynomial.legendre.leggauss(np.shape(old_psi)[1])[1]
    if(ext_source!=None):    
        for i in range(np.size(source)):
            source[i]=0.5*(Sigma_s*np.dot(weights,old_psi[i])+ext_source)     #source = 1/2(Sigma_s*flux+Q), weighted sum is a scalar product
    else:
        for i in range(np.size(source)):
            source[i]=0.5*(Sigma_s*np.dot(weights,old_psi[i])) 
    #new flux calculation
    points=np.polynomial.legendre.leggauss(np.shape(old_psi)[1])[0]
    new_psi=np.zeros(np.shape(old_psi))
    for i in range(np.shape(old_psi)[1]):        #number of columns (discrete ordinates)
        if(points[i]>=0):
            if(zeroth==True):
                psi_minus=boundary_conditions[0](points[i])
            else:
                psi_minus=0
            for j in range(np.shape(old_psi)[0]):
                new_psi[j][i]=(source[j]+psi_minus*(2*points[i]*J/a))/(2*points[i]*J/a + Sigma_t)   #calculate new angular flux value
                psi_minus=2*new_psi[j][i]-psi_minus     #get new boundary value (psi_plus) for next spatial segment
        else:
            if(zeroth==True):            
                psi_plus=boundary_conditions[1](points[i])
            else:
                psi_plus=0
            for j in reversed(range(np.shape(old_psi)[0])):
                new_psi[j][i]=(source[j]-psi_plus*(2*points[i]*J/a))/(-2*points[i]*J/a + Sigma_t)   #the denominator is always positive because the cosine is negative
                psi_plus=2*new_psi[j][i]-psi_plus     #get new boundary value (psi_minus) for next spatial segment
        
    return new_psi



    
    



