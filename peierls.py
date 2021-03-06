# Various functions for handling Peierls potentials 

from numpy import *
from scipy import interpolate, integrate
import peierls_pot as pot

def num_integ(inx,iny):
    """Given values of a function (iny) evaluated at points (inx), 
       calculate integral of function at points using trapezium rule"""
    inx = atleast_1d(inx)
    out = zeros(inx.shape, dtype=inx.dtype)
    for n in xrange(len(out)):
        out[n] = integrate.trapz(iny[0:n],inx[0:n])
    return out


def get_u0_element (disp, potderiv, stress):
    u_0 = 0.0
    u_0_index = 0
    for n in xrange(len(disp)):
        if (potderiv[n] > stress):
            u_0 = disp[n]
            u_0_index = n
            break

    return(u_0, u_0_index)

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
        (u_0[i], u_0_index[i]) = get_u0_element (disp, potderiv, stresses[i])

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

def get_u_max_from_h(x, h, u_0, u_0_index):
    """Given a Peiels potential x, a kink height
    h, and a initial displacment u_0 calculate
    the index of the kink heigh u_max_index and 
    the actual position of the kink u_max. For
    e.q. 3 of Carrez 2009"""

    u_max = max(x)
    u_max_index = len(x)
    for i in xrange(len(x)):
        if (x[i] < u_0):
            continue
        if (x[i] > u_0+h):
            u_max = u_0 + h
            u_max_index = i
            break

    return(u_max, u_max_index)

def square_dislo_energy_screw_elas_element(x, y, h, w, roh, shear_mod, poss, burgers):
    # screw dislocation with edge kinks. 
    # EQ. 4 of Carrez et al. 09
    Eelas = (shear_mod * burgers**2) / ( 2.0 * pi) * \
            ( sqrt(w**2+h**2)-w-h \
              + (w*log((2.0*w)/(w+sqrt(w**2+h**2)))) \
              - ((1.0/(1.0-poss))*(w-sqrt(w**2+h**2)+h*log((h+sqrt(w**2+h**2))/w) - \
              h*log(h/(roh*e)))) ) # NB - e is 2.71... from numpy. (not exp(roh) !!)
    return Eelas

def square_dislo_energy_peierls_element(x, y, yderv, h, w, stress):
    (u_0, u_0_index) = get_u0_element(x, yderv, stress)
    (u_max, u_max_index) = get_u_max_from_h(x, h, u_0, u_0_index)
    if (h+u_0 < max(x)):
        Epeierls = 2 * num_integ(x[u_0_index+1:u_max_index], y[u_0_index+1:u_max_index])[-1] + \
                     w * (y[u_max_index] - y[u_0_index])
    else:
        Epeierls = None
    return Epeierls

def square_dislo_energy_work_element(stress, burgers, h, w):
    return (stress * burgers * h * w)

 
def square_dislo_energy_screw_element(x, y, yderv, h, w, roh, stress, shear_mod, poss, burgers):
    Epeierls = square_dislo_energy_peierls_element(x, y, yderv, h, w, stress)
    if Epeierls is None:
        Epeierls = 0.0
        Eelas = 0.0
        Ework = 0.0
        Etot = 0.0
    else:
        Eelas = square_dislo_energy_screw_elas_element(x, y, h, w, roh, shear_mod, poss, burgers)
        Ework = square_dislo_energy_work_element(stress,burgers,h,w)
        Etot = Eelas + Epeierls - Ework
    return (Etot, Eelas, Epeierls, Ework)

def square_dislo_energy_screw_deriv(x, y, yderv, h, w, roh, stress, shear_mod, poss, burgers, dh, dw):
    (Etot, Eelas, Epeierls, Ework) = square_dislo_energy_screw_element(x, y, yderv, h, w, roh, stress, shear_mod, poss, burgers)
    (Eph, a, b, c) = square_dislo_energy_screw_element(x, y, yderv, h+dh, w, roh, stress, shear_mod, poss, burgers)
    (Emh, a, b, c) = square_dislo_energy_screw_element(x, y, yderv, h-dh, w, roh, stress, shear_mod, poss, burgers)
    (Epw, a, b, c) = square_dislo_energy_screw_element(x, y, yderv, h, w+dw, roh, stress, shear_mod, poss, burgers)
    (Emw, a, b, c) = square_dislo_energy_screw_element(x, y, yderv, h, w-dw, roh, stress, shear_mod, poss, burgers)
    dEbydH = (Eph - Emh) / 2 * dh
    dEbydW = (Epw - Emw) / 2 * dw
    return (Etot, Eelas, Epeierls, Ework, dEbydH, dEbydW)

    

def square_dislo_energy_screw (x, y, yderv, h, w, roh, stress, shear_mod, poss, burgers, dH=None, dW=None):

    Epeierls = zeros((len(stress),len(h),len(w)), dtype=stress.dtype)
    Eelas = zeros((len(stress),len(h),len(w)), dtype=stress.dtype)
    Ework = zeros((len(stress),len(h),len(w)), dtype=stress.dtype)
    Etot = zeros((len(stress),len(h),len(w)), dtype=stress.dtype)
    dEbydW = zeros((len(stress),len(h),len(w)), dtype=stress.dtype)
    dEbydH = zeros((len(stress),len(h),len(w)), dtype=stress.dtype)
    if ((dH is None) and (dW is None)):
        for i in xrange(len(stress)):
            for j in xrange(len(h)):
                for k in xrange(len(w)):
                    (Etot[i,j,k], Eelas[i,j,k], Epeierls[i,j,k], Ework[i,j,k]) \
                       = square_dislo_energy_screw_element(x, y, yderv, h[j], w[k], roh, stress[i], shear_mod, poss, burgers)
    elif ((dH is None) or (dW is None)):
        raise ValueError("Cannot handle deriv subsets")
    else:
        for i in xrange(len(stress)):
            for j in xrange(len(h)):
                for k in xrange(len(w)):
                    (Etot[i,j,k], Eelas[i,j,k], Epeierls[i,j,k], Ework[i,j,k], dEbydH[i,j,k], dEbydW[i,j,k]) \
                       = square_dislo_energy_screw_deriv(x, y, yderv, h[j], w[k], roh, stress[i], shear_mod, poss, burgers, dH, dW)
        
    
    # Answer in eV
    Na = 6.022E23 # Avagadro's number
    Etot = ((Etot * Na)/1000) / 96
    Eelas = ((Eelas * Na)/1000) / 96
    Epeierls = ((Epeierls * Na) / 1000) / 96
    Ework = ((Ework * Na) / 1000) / 96
    return (Etot, Eelas, Epeierls, Ework, dEbydH, dEbydW)

    


def kinker (x, y, G=None, b=None, silent=False, method=None, params=None, maxx=None):
    """Solve the line tension model for a set of points 
       describing the Peiels potential. Input is:
       x: list of displacments from 0 to b
       y: list of energies.
       G: optional bulk modulous to estimate elastic
          contribution to energy.
       b: optional value of Burgers vector (otherwise assume max(x)
       Output is:
       xfine: interpolation over x
       yfine: interpolation over y
       sigma_p: The Peiels stress
       zdiff: shape of geometrical kink (one point per xfine)
       """

    Na = 6.022E23 # Avagadro's number

    if maxx is None:
        #we need to know the maximum position of the x axis, which we take as the 
        # periodicity of the potential 
        x_max = max(x)
        if not silent:
            print "From the data it looks like x max' is %5g m" % x_max
        if b is None:
            b = x_max
    else:
        x_max = maxx

    # Fit to cubic spline - avoid having xfine = 0 as this gives 
    # NaN when doing 1/sqrt(y) also because close to the origin we 
    # can have yfine[1] < yfine[0] (for example)
    xfine = arange(0.025E-10,(x_max - 0.025E-10),0.001E-10)
    if method is None:
        tck = interpolate.splrep(x,y,s=0)
        yfine = interpolate.splev(xfine,tck,der=0)
        # Calculate derivative 
        yderfine = interpolate.splev(xfine,tck,der=1)
    else:
        # Probably should pass in the func, which would be better and could
        # generalize to interp too?
        yfine = pot.eval_func(xfine,params,x_max,method)
        # Calculate derivative 
        yderfine = pot.pot_num_deriv(xfine,params,x_max,(0.025E-10/10),method)

    # For PN data we'll need to add the elastic energy of the dislocation.
    # Approximatly: 1/2Gb^2
    # FIXME: caller should do this!
    if G != None:
        yfine = yfine + (0.5 * G * x_max**2.0)

    # Calculate Sd...
    Sd = average(yfine) 
    if not silent:
        print "Sd has a value of %10.10g Jm-1" % Sd

    if not silent:
        print "Maximum peierls pot is %5g Jm-1" % max(yfine)
        print "Maximum peierls pot is %5g 10E-3 eV.Ang-1" % (max(yfine)*6.24150974E8/10E-3)
        print "Minimum peierls pot is %5g Jm-1" % min(yfine)
        print "Minimum peierls pot is %5g 10E-3 eV.Ang-1" % (min(yfine)*6.24150974E8/10E-3)
        print "Vp (max - min) is %5g Jm-1" % (max(yfine)-min(yfine))
        print "Vp (max - min) is %5g 10E-3 eV.Ang-1" % ((max(yfine)-min(yfine))*6.24150974E8/10E-3)


    sigma_p_index = argmax(yderfine)
    sigma_p = yderfine[sigma_p_index] / x_max
    if not silent:
        print "Sigma_p has a value of %5g Pa " % sigma_p

    w_k = x_max * sqrt(Sd/(2.0*(yfine[argmax(yfine)]-yfine[0])))

    if not silent:
        print "w_k has a calculated value of %5g m " % w_k

    # Calculate 1/sqrt(U) and plot
    invsqrt_yfine = 1.0/(sqrt(yfine-yfine[0]))
    invsqrt_y = 1.0/(sqrt(y-y[0]))

    # Calculate single kink shape - equation 5 
    # of Seeger 1981. Note that U(0) = 0. (no applied stress)

    zdiff = zeros(xfine.shape, dtype=xfine.dtype)
    zdiff = sqrt(Sd/2.0) * num_integ(xfine[1:],(1.0/(sqrt(yfine[1:]-yfine[0])))) 
    kink_energy = sqrt(2.0*Sd) * (num_integ(xfine,sqrt(yfine-yfine[0])))[-1] 
    if not silent:
        print "Energy of a geometrical kink = %10.10g J " % kink_energy 
        print "Energy of a geometrical kink = %10.10g kJ/mol " % ((kink_energy * Na)/1000)

    # Calculate kink energy according to Dorn & Rajnak 
    kink_energy2 = yfine[0] * (num_integ(xfine,sqrt(((yfine/yfine[0])**2) - 1  )))[-1] 
    if not silent:
        print "Energy of a geometrical kink (2) = %10.10g J " % kink_energy2
        print "Energy of a geometrical kink = %10.10g kJ/mol " % ((kink_energy2 * Na)/1000)

    sigma = array([0.001*sigma_p, 0.01*sigma_p, 0.1*sigma_p, 0.2*sigma_p, 0.25*sigma_p, \
              0.3*sigma_p, 0.4*sigma_p, \
              0.5*sigma_p, 0.6*sigma_p, 0.7*sigma_p, \
              0.75*sigma_p, 0.8*sigma_p,  0.9*sigma_p, 0.99*sigma_p, 0.999*sigma_p])

    # kp energy takes sigma * burgers vector as input (c.f. seeger e.q. 9) but I think 
    # this is actually the periodicty of the function and Seeger assumes a edge dislo.
    (u_0, u_max, H_kp, Un, zdiff_kp) = kp_energy(xfine, yfine, yderfine, (sigma*x_max))


    if not silent:
        print     "Sigma (Pa)     u_0 (m)        u_max (m)   H_kp (kJ/mol) Un"
        for i in xrange(len(sigma)):
            print "% .3g    % .3g      % .3g    % .3g    % .3g" % (sigma[i], u_0[i], u_max[i], (H_kp[i]*Na/1000), (Un[i]*Na/1000))

    return (xfine, yfine, sigma_p, zdiff, u_0, u_max, H_kp, Un, zdiff_kp, (yderfine/x_max), sigma, sigma_p_index, kink_energy2)
