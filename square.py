#!/usr/bin/python

from peierls import square_dislo_energy_screw
import numpy as np
import scipy as sp
import matplotlib 
from pylab import *

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
    x,y = np.loadtxt(FILE, unpack=True)
    FILE.close()
    return(x,y)


basename = raw_input("Basename (input is basname.dat): ")
if (basename == ""):
    basename = 'example'
    filename = None
else:
    filename = basename + '.dat'

xmax_input = raw_input("Provide a maximum value of x in m (empty string to use input file and 5.0E-10 for old results):")
if (xmax_input == ""):
    xmax_input = None
else:
    xmax_input = float(xmax_input)


(x,y) = load_data(filename)

if xmax_input is None:
    xmax = max(x)
else:
    xmax = xmax_input

xfine = np.arange(0.025E-10,(xmax - 0.025E-10),0.001E-10)
tck = sp.interpolate.splrep(x,y,s=0)
yfine = sp.interpolate.splev(xfine,tck,der=0)
yderv = sp.interpolate.splev(xfine,tck,der=1)

burgers = sqrt(xmax**2+xmax**2)
h = np.array([0.12, 0.17, 0.22, 0.27, 0.32, 0.37, 0.42, 0.47, 0.52, 0.57, 0.62, 0.67, 0.72, 0.77, 0.82, 0.87, 0.92])
#h = np.array([0.1, 0.52])
h = h * xmax
w = np.array([200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 2000, 3000, 4000])
w = w * burgers
roh = 0.05 * 2.08 * xmax 
stress = np.array([1.0E5, 10E6, 15E6, 20E6, 25E6, 30E6])
poss = 0.18
shear_mod = 116.5E9


(Etot, Eelas, Epeierls, Ework) = square_dislo_energy_screw (xfine, yfine, yderv, h, w, roh, stress, shear_mod, poss, burgers)

print "Etot:"
print Etot
print "Eelas:"
print Eelas
print "Ep:"
print Epeierls
print "Work:"
print Ework

# plot the figures
for i in xrange(len(stress)):
    figure(i)
    title(str(stress[i])+ " Pa")
    contour(w,h,Etot[i,:,:], 40)
    colorbar()
show()

