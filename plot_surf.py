#!/usr/bin/python

from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate
from peierls import square_dislo_energy_screw

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


basename = raw_input("Basename (input is basname.dat): ")
if (basename == ""):
    basename = 'example'
    filename = None
else:
    filename = basename + '.dat'


# Load data set, x is the displacment (u) and y is the energy (U)
#(x,y) = load_data(None) # For test data
(x,y) = load_data(filename)

G = 100.0E9
y = y + (0.5 * G * 4.498641E-20)

Na = 6.022E23 # Avagadro's number

# Fit to cubic spline - avoid having xfine = 0 as this gives 
# NaN when doing 1/sqrt(y) also because close to the origin we 
# can have yfine[1] < yfine[0] (for example)
#xfine = arange(0.025E-10,4.975E-10,0.001E-10)
xfine = arange(0.025E-10,2.000E-10,0.001E-10)
tck = interpolate.splrep(x,y,s=0)
yfine = interpolate.splev(xfine,tck,der=0)

# Plot data 
plt.figure(1)
plt.plot(x,y,'o',xfine,yfine,'-')
plt.ylabel('U')
plt.xlabel('u')
plt.title('U(u)')

# Calculate derivative 
yderfine = interpolate.splev(xfine,tck,der=1)
sigma_p_index = argmax(yderfine)
sigma_p = yderfine[sigma_p_index]
#print "Sigma_p has a value of %5g Pa " % (sigma_p / 5.0E-10)
print "Sigma_p has a value of %5g Pa " % (sigma_p / 2.121E-10)
print "Sigma_p has a value of %5g MPa " % (sigma_p*1.0E-6 / 2.08E-10)


delta_h = 0.025
delta_w = 2.500

hs = arange(0, 1, delta_h)
ws = arange(0, 5000, delta_w)

Etot = square_dislo_energy_screw(xfine, yfine, hs[:], ws[:], 0.2, 0.01,10,10,10,10)
