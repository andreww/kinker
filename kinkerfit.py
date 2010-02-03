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

    calc_x = (sigma_b/sigma_p)
    calc_y = (Un/(kink_energy2*2.0))
    interp_calc_y = interp(T_n,calc_x,calc_y)
    error = (tau_n - interp_calc_y)
    print "----------------------------------------"
    print "kink energy" 
    print  kink_energy2
    print "p stress" 
    print  sigma_p
    print "sum squares"
    print  sum((error*error))
    print "time" 
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
p0 = [2.4E-10, 2.4E-10, 1, 2, 1, 2, 3, 4, 1, 2]
opt_p = pot.fit_pot(p0, x, y, x_len)


# Initial solution with this function
(xfine, yfine, sigma_p, zdiff, u_0, u_max, H_kp, Un, \
 zdiff_kp, yderfine, sigma_b, sigma_p_index, kink_energy2) = kinker(x,y,G=60E9,method='func',params=opt_p)

# Load experimental data
exptdata = raw_input("Expt data (input is basname.dat): ")
(expt_T,expt_tau) = pot.load_data(exptdata)

T_n = expt_T/500.0
# Don't know why I need to / 5000
# but 9.0 is tau at Tcrit (removes background
# from this dataset).
tau_n = (expt_tau-9.0)/5000

print "Optimizing potential..."
print "Kink energy, sigma_p, sumsq, time"
t0 = time.time()
#full_opt_p, sucess = optimize.leastsq(errfunc, opt_p, args=(x, y))

full_opt_p = [  2.40000000e-10,   2.40000000e-10,  -1.15171225e-03,   1.69230101e-03, \
                     1.38187764e-06,   2.00000000e+00,   3.00000000e+00,   5.70605907e-02, \
                     8.88014468e+00,   2.03026690e+00]
sucess = 10

print "Optimal parameter set and sucess variable:"
print full_opt_p
print sucess

print "Solution to kink model with this paramerer set"
(xfine, yfine, sigma_p, zdiff, u_0, u_max, H_kp, Un, \
 zdiff_kp, yderfine, sigma_b, sigma_p_index, kink_energy2) = kinker(x,y,G=60E9,method='func',params=full_opt_p)

# Plot data 
plt.figure(1)
plt.plot(x,y,'o',x,pot.pot_func(x,opt_p, x_len),'-',x,pot.pot_func(x,full_opt_p, x_len),'--')
plt.legend(('Starting points', 'Fit to starting points', 'fit to expt data'))
plt.ylabel('U')
plt.xlabel('u')
plt.title('U(u)')

plt.figure(2)
plt.plot(xfine[1:],zdiff,'-')
plt.ylabel('z-z0')
plt.xlabel('u')
plt.title('Geometrical kink shape')

plt.figure(3)
plt.plot(xfine,yderfine,'--',u_0,sigma_b,'o',u_max,sigma_b,'+',array([xfine[sigma_p_index]]),array([sigma_p]),'x')
plt.title('Derivative of potential')
plt.ylabel('dU/du')
plt.xlabel('u')

plt.figure(4)
plt.plot((sigma_b/5.0E-10),(Un*Na/1000),'--')
plt.title('Critical stress / kink energy')
plt.xlabel('Tau* - Pa')
plt.ylabel('Un - kJ/mol')

plt.figure(5)
plt.plot((sigma_b/sigma_p),(Un/(kink_energy2*2.0)),'--',(tau_n/sigma_p),T_n,'o')
plt.title('Critical stress / kink energy')
plt.xlabel('Tau*/sigma_p')
plt.ylabel('Un/2Uk')

plt.show()
