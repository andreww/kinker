# Various functions for handling Peierls potentials 

from numpy import *
from scipy import interpolate, integrate

def num_integ(inx,iny):
    """Given values of a function (iny) evaluated at points (inx), 
       calculate integral of function at points using trapezium rule"""
    inx = atleast_1d(inx)
    out = zeros(inx.shape, dtype=inx.dtype)
    for n in xrange(len(out)):
        out[n] = integrate.trapz(iny[0:n],inx[0:n])
    return out


def get_u0 (disp, potderiv, stresses):
    """Given the derivative of the Peierls potential 
       as two arrays (dispm, potderiv) and an array 
       of applied stresses calculate the position of 
       the straight dislocation on the potential surface.
       i.e. solve equation 4 of Seeger or 2 of Dorn and 
       Rajnak"""
    u_0 = zeros(stresses.shape, dtype=stresses.dtype)
    u_0_index = zeros(stresses.shape, dtype=integer) # What index do we want to start from?
    for i in xrange(len(stresses)):
        for n in xrange(len(disp)):
            if (potderiv[n] > stresses[i]):
                u_0[i] = disp[n]
                u_0_index[i] = n
                break

    return(u_0, u_0_index)

def kp_energy (xfine, yfine, yderfine, sigma_b):
    """Calculate the kink pair formation energy from Seegers model 
       using the equations of Seeger or Dorn and Rajnak (they differ)
       from the Peierls potental (xfine, yfine) and it's derivative
       at a set of stresses (sigma_b).
       Returns a tuple of (u_0, u_max, H_kp, Un)"""

    # Calculate Sd...
    Sd = average(yfine)

    (u_0, u_0_index) = get_u0(xfine, yderfine, sigma_b)

    # Upper limits (u_max) come from eq.9: this is a brute force search.
    
    u_max = zeros(sigma_b.shape, dtype=sigma_b.dtype)
    u_max_index = zeros(sigma_b.shape, dtype=integer)
    U_u_max = zeros(sigma_b.shape, dtype=sigma_b.dtype)
    inner_func = zeros((len(sigma_b),len(xfine)))
    inner_func2 = zeros((len(sigma_b),len(xfine)))
    
    for i in xrange(len(sigma_b)):
        # Build an array inner_func = U(u) - U(u0) - (u - u_0)b_sigma
        # First index over values of b_sigma and second over values of u.
        # This is the inner bit of equation 8 and the LHS of equation 9.
        inner_func[i,:] = (yfine - yfine[u_0_index[i]] - ((xfine-u_0[i])*sigma_b[i]))
        inner_func2[i,:] = (yfine**2 - (((xfine-u_0[i])*sigma_b[i]) + yfine[u_0_index[i]])**2)
    
        for n in xrange(len(xfine)):
            if (n < (u_0_index[i]+5)):
                continue # u_max must be after u_0 and sometimes these seems to be a numerical
                         # issue so make it at least 5 spaces after u_0
            if (inner_func[i,n] < 0.0):
                u_max[i] = xfine[n]
                u_max_index[i] = n
                U_u_max[i] = yfine[n]
                break
    
    # Calculate kink pair formation v's stress
    H_kp = zeros(sigma_b.shape, dtype=sigma_b.dtype)
    Un = zeros(sigma_b.shape, dtype=sigma_b.dtype)
    for i in xrange(len(sigma_b)):
        H_kp[i] = 2.0 * sqrt(2.0*Sd) * \
                    (num_integ(    \
                      xfine[(u_0_index[i]):u_max_index[i]], \
                      sqrt(inner_func[i,(u_0_index[i]):u_max_index[i]]))[-1])
    
        Un[i] = 2.0 * (num_integ(xfine[(u_0_index[i]):u_max_index[i]], \
                      sqrt(inner_func2[i,(u_0_index[i]):u_max_index[i]]))[-1])

    # Calculate shape of kink pair - this looks odd? 
    
    zdiff_kp = zeros((len(sigma_b),len(xfine)))
    for i in xrange(len(sigma_b)):
        zdiff_kp[i,(u_0_index[i]+1):u_max_index[i]] = (sqrt(Sd/2.0) * \
                     (num_integ(xfine[(u_0_index[i]+1):u_max_index[i]], \
                     (1.0/(sqrt(yfine[(u_0_index[i]+1):u_max_index[i]]-yfine[u_0_index[i]]))))))


    return(u_0, u_max, H_kp, Un, zdiff_kp)

