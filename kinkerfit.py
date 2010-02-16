#!/usr/bin/python

from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate, optimize
from peierls import kp_energy, num_integ, kinker
import time as time
import peierls_pot as pot

def errfunc(p, x, y):
    # Assume we have loads of cal data so we can
    # intepolate it to find values equiv to the 
    # experimental points (using numpys numpy.interp)
    (xfine, yfine, sigma_p, zdiff, u_0, u_max, H_kp, Un, \
     zdiff_kp, yderfine, sigma_b, sigma_p_index, kink_energy2) = kinker(x,y, G=60E9,silent=True,method='func',params=p )

    calc_y = (Un/(kink_energy2*2.0))*500.0
    calc_x = sigma_b 
    interp_calc_y = interp(tau_n,calc_x,calc_y)
    error = (T_n - interp_calc_y)
    print "----------------------------------------"
    print "kink energy (J)" 
    print  kink_energy2
    print "p stress (Pa)" 
    print  sigma_p
    print "sum squares (K^2)"
    print  sum((error*error))
    print "time (s)" 
    print (time.time() - t0)
    print "----------------------------------------"
    return error

Na = 6.022E23 # Avagadro's number

# Load data set, x is the displacment (u) and y is the energy (U)
#(x,y) = load_data(None) # For test data
basename = raw_input("Basename (input is basname.dat): ")
if (basename == ""):
    basename = 'example'
    filename = None
else:
    filename = basename + '.dat'

(x,y) = pot.load_data(filename)
x_len = max(x)

# Solve for inital potential using bspline interpolation...
print "Solving kink nucleation problem with inital potential using bspline interp"
(xfine_interp, yfine_interp, sigma_p, zdiff, u_0, u_max, H_kp, Un, \
 zdiff_kp, yderfine, sigma_b, sigma_p_index, kink_energy2) = kinker(x,y,G=60E9)

# Now do it again using a sine function...

# Fit this data to a sin function
print "Fitting initial potential to sin function"
#p0 = [2.4E-10, 2.4E-10, 1, 2, 1, 2, 3, 4, 1, 2]
p0 = [2.4, 0.001, 0.001, 2.0, 2.0, 2.0]
opt_p = pot.fit_pot(p0, x, y, x_len)


# Initial solution with this function
(xfine_init, yfine, sigma_p, zdiff_init, u_0_init, u_max, H_kp, Un_init, \
 zdiff_kp, yderfine_init, sigma_b_init, sigma_p_index, kink_energy2_init) = kinker(x,y,G=60E9,method='func',params=opt_p)

# Load experimental data
exptdata = raw_input("Expt data (input is basname.dat): ")
(T_n,expt_tau) = pot.load_data(exptdata)
# Convert from MPa to Pa and subtract tau_crit (15 MPa, T_crit is 500 K, handled in error func.)
tau_n = (expt_tau-15.0)*1E6

print "Optimizing potential..."
t0 = time.time()
full_opt_p, cov_x, infodict, mesg, sucess = optimize.leastsq(errfunc, opt_p, args=(x, y), full_output=1)
print mesg
print "Optimal parameter set:"
print full_opt_p
print "sucess code:"
print sucess

print "Solution to kink model with this paramerer set"
(xfine, yfine, sigma_p, zdiff, u_0, u_max, H_kp, Un, \
 zdiff_kp, yderfine, sigma_b, sigma_p_index, kink_energy2) = kinker(x,y,G=60E9,method='func',params=full_opt_p)

# Plot data 
plt.figure(1)
plt.plot(x,y,'o',x,pot.pot_func(x,opt_p, x_len),'-',x,pot.pot_func(x,full_opt_p, x_len),'--')
plt.legend(('Starting points', 'Fit to starting points', 'fit to expt data'))
plt.ylabel('U (J/m)')
plt.xlabel('u (m)')
plt.title('U(u)')

plt.figure(2)
plt.plot(xfine[1:],zdiff,'-',xfine_init[1:],zdiff_init,'--')
plt.ylabel('z-z0 (m)')
plt.xlabel('u (m)')
plt.legend(('fit to expt data', 'init model'))
plt.title('Geometrical kink shape')

plt.figure(3)
plt.plot(xfine,yderfine,'-',u_0,sigma_b,'o',xfine_init,yderfine_init,'--',u_0_init,sigma_b_init,'x')
plt.title('Derivative of potential')
plt.ylabel('dU/du (MPa)')
plt.xlabel('u (m)')
plt.legend(('fit to expt data', 'stresses used for expt data', 'init model', 'init model stresses'))

plt.figure(4)
plt.plot((sigma_b/5.0E-10),(Un*Na/1000),'--')
plt.title('Critical stress / kink energy')
plt.xlabel('Tau* - Pa')
plt.ylabel('Un - kJ/mol')

plt.figure(5)
plt.plot(sigma_b,((Un/(kink_energy2*2.0))*500),'--',(tau_n),T_n,'o',sigma_b_init,((Un_init/(kink_energy2_init*2.0))*500),'-')
plt.legend(('fit to expt data', 'expt data', 'init model'))
plt.title('Critical stress / kink energy')
plt.xlabel('Tau* (MPa)')
plt.ylabel('Un/2Uk*T_crit (K)')

plt.show()
