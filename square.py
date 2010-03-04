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
h = np.array([0.02, 0.07, 0.12, 0.17, 0.22, 0.27, 0.32, 0.37, 0.42, 0.47, 0.52, 0.57,  0.62, 0.67, 0.72, 0.77, 0.82, 0.87, 0.92, 0.97])
#h = np.array([0.1, 0.52])
h = h * xmax
w = np.array([50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500])
w = w * burgers
roh = 0.05 * 2.08 * xmax 
#stress = np.array([1.0E5, 10E6, 15E6, 20E6, 25E6, 30E6])
stress = np.array([20E6])
poss = 0.18
shear_mod = 116.5E9


(Etot, Eelas, Epeierls, Ework, dEbydH, dEbydW) = square_dislo_energy_screw (xfine, yfine, yderv, h, w, roh, stress, shear_mod, poss, burgers, dH=0.01E-10, dW=1.0E-10)

deriv = abs(dEbydH) + abs(dEbydW)
#derivderiv = np.zeros((len(stress),10,10), dtype=stress.dtype)
#for i in xrange(len(stress)):
#    deriv_func = sp.interpolate.interp2d(w, h, deriv[i,:,:])
#
#
#    hh = np.linspace(min(h), max(h), 10)
#    ww = np.linspace(min(w), max(w), 10)
#    derivderiv[i,:,:] = deriv_func(ww, hh)
#

print "Etot:"
print Etot
print "Eelas:"
print Eelas
print "Ep:"
print Epeierls
print "Work:"
print Ework
print "dEbydH:"
print dEbydH
print "dEbydW:"
print dEbydW
print "|dE/dH| + |dE/dW|:"
print deriv

allH = np.zeros(deriv.shape)
allW = np.zeros(deriv.shape)
print allW.shape
print deriv.shape
print w.shape
for i in xrange(len(stress)):
    for j in xrange(len(w)):
        for k in xrange(len(h)):
            # Don't really understand the inversion 
            # of axes here...
            allW[i,j,k] = w[k]
            allH[i,j,k] = h[j]

#print deriv.argmin(axis=0)
Ecrit = [0.00]*len(stress)
Hcrit = [0.00]*len(stress)
Wcrit = [0.00]*len(stress)
for i in xrange(len(stress)):
    lin_deriv = deriv[i,:,:].ravel()
    lin_allW = allW[i,:,:].ravel()
    lin_allH = allH[i,:,:].ravel()
    lin_Etot = Etot[i,:,:].ravel()
    pos = lin_deriv.argmin()
    Ecrit[i] = lin_Etot[pos]
    Wcrit[i] = lin_allW[pos]/1E-10
    Hcrit[i] = lin_allH[pos]/1E-10
    print 'Stress: ' + str(stress[i])
    print 'Ecrit: ' + str(Ecrit[i])
    print 'Wcrit: ' + str(Wcrit[i])
    print 'Hcrit: ' + str(Hcrit[i])


# plot the figures
fnum = 0
for i in xrange(len(stress)):
    fnum = fnum + 1
    figure(fnum)
    title(str(stress[i]*1E-6)+ " MPa")
    contour(w/1E-10,h/1E-10,Etot[i,:,:], 40)
    colorbar()
    scatter(Wcrit[i], Hcrit[i], marker='o')
    scatter(allW[i,:,:]/1E-10, allH[i,:,:]/1E-10, marker='x')
    fnum = fnum + 1
    figure(fnum)
    title("derivative of " + str(stress[i]*1E-6)+ " Pa")
    contour(w/1E-10,h/1E-10,deriv[i,:,:], 20)
    colorbar()
    scatter(Wcrit[i], Hcrit[i], marker='o')
    scatter(allW[i,:,:]/1E-10, allH[i,:,:]/1E-10, marker='x')
show()

