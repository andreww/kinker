#!/usr/bin/python

from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate
from peierls import kp_energy, num_integ, kinker

def load_data (filename):
    """Load data from filename, if None load example data."""
    from StringIO import StringIO   # StringIO behaves like a file object
    example_data = "0.0001E-10 1.0E-9\n0.5E-10 1.2E-9\n1.0E-10 1.5E-9\n1.5E-10 1.6E-9\n2.0E-10 1.7E-9\n \
                    2.5E-10 1.8E-9\n3.0E-10 1.7E-9\n3.5E-10 1.6E-9\n4.0E-10 1.5E-9\n4.5E-10 1.2E-9\n4.9999E-10 1.0E-9"
    if (filename==None):
        FILE = StringIO(example_data)
        print "Using example data"
    else:
        FILE = open(filename,'r')
        print "Loading data from %s" % filename
    x,y = loadtxt(FILE, unpack=True)
    FILE.close()
    return(x,y)


def errorfunc(p, x, y):
    # Assume we have loads of cal data so we can
    # intepolate it to find values equiv to the 
    # experimental points (using numpys numpy.interp)
    ## FIXME - I'm stuck at the mo
    (xfine, yfine, sigma_p, zdiff, u_0, u_max, H_kp, Un, \
     zdiff_kp, yderfine, sigma_b, sigma_p_index, kink_energy2) = kinker(x,y)
    interp_calc_y = interp(x, , calc_y)
    error = (expt_y - interp_calc_y)
    return error

basename = raw_input("Basename (input is basname.dat): ")
if (basename == ""):
    basename = 'example'
    filename = None
else:
    filename = basename + '.dat'

Na = 6.022E23 # Avagadro's number

# Load data set, x is the displacment (u) and y is the energy (U)
#(x,y) = load_data(None) # For test data
(x,y) = load_data(filename)

# Initial solution with out opt.
(xfine, yfine, sigma_p, zdiff, u_0, u_max, H_kp, Un, \
 zdiff_kp, yderfine, sigma_b, sigma_p_index, kink_energy2) = kinker(x,y)

FILE = open('MgO_experiments.dat', 'r')
expt_T,expt_tau = loadtxt(FILE, unpack=True)
FILE.close()

T_n = T/500.0
# Don't know why I need to / 5000
# but 9.0 is tau at Tcrit (removes background
# from this dataset).
tau_n = (tau-9.0)/5000

opt_y, sucess = optimize.leastsq(errfunc, x, args=(x, y))




# Plot data 
plt.figure(1)
plt.plot(x,y,'o',xfine,yfine,'-')
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
plt.plot((sigma_b/sigma_p),(Un/(kink_energy2*2.0)),'--',(sigma_b/sigma_p),(Un/(kink_energy2*2.0)),'o')
plt.title('Critical stress / kink energy')
plt.xlabel('Tau*/sigma_p')
plt.ylabel('Un/2Uk')

plt.show()
