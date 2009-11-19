#!/usr/bin/python

from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate
from peierls import kp_energy, num_integ

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

Na = 6.022E23 # Avagadro's number

# Fit to cubic spline - avoid having xfine = 0 as this gives 
# NaN when doing 1/sqrt(y) also because close to the origin we 
# can have yfine[1] < yfine[0] (for example)
xfine = arange(0.025E-10,4.975E-10,0.001E-10)
tck = interpolate.splrep(x,y,s=0)
yfine = interpolate.splev(xfine,tck,der=0)

# Plot data 
plt.figure(1)
plt.plot(x,y,'o',xfine,yfine,'-')
plt.ylabel('U')
plt.xlabel('u')
plt.title('U(u)')

# Calculate Sd...
Sd = average(yfine) 
print "Sd has a value of %10.10g Jm-1" % Sd

# Calculate derivative 
yderfine = interpolate.splev(xfine,tck,der=1)
sigma_p_index = argmax(yderfine)
sigma_p = yderfine[sigma_p_index]
print "Sigma_p has a value of %5g Pa " % (sigma_p / 5.0E-10)

w_k = 5.0E-10 * sqrt(Sd/(2.0*(yfine[argmax(yfine)]-yfine[0])))

print "w_k has a calculated value of %5g m " % w_k

# Calculate 1/sqrt(U) and plot
invsqrt_yfine = 1.0/(sqrt(yfine-yfine[0]))
invsqrt_y = 1.0/(sqrt(y-y[0]))

plt.figure(2)
# NB - avoid plotting end points (infnty)
plt.plot(x[1:-1],invsqrt_y[1:-1],'o',xfine[1:-1],invsqrt_yfine[1:-1],'-')
plt.ylabel('U^-1/2')
plt.xlabel('u')
plt.title('U(u)^-1/2')

# Calculate single kink shape - equation 5 
# of Seeger 1981. Note that U(0) = 0. (no applied stress)

zdiff = zeros(xfine.shape, dtype=xfine.dtype)
zdiff = sqrt(Sd/2.0) * num_integ(xfine[1:],(1.0/(sqrt(yfine[1:]-yfine[0])))) 
kink_energy = sqrt(2.0*Sd) * (num_integ(xfine,sqrt(yfine-yfine[0])))[-1] 
print "Energy of a geometrical kink = %10.10g J " % kink_energy 
print "Energy of a geometrical kink = %10.10g kJ/mol " % ((kink_energy * Na)/1000)

# Calculate kink energy according to Dorn & Rajnak 
kink_energy2 = yfine[0] * (num_integ(xfine,sqrt(((yfine/yfine[0])**2) - 1  )))[-1] 
print "Energy of a geometrical kink (2) = %10.10g J " % kink_energy2
print "Energy of a geometrical kink = %10.10g kJ/mol " % ((kink_energy2 * Na)/1000)

plt.figure(3)
plt.plot(xfine[1:],zdiff,'-')
plt.ylabel('z-z0')
plt.xlabel('u')
plt.title('Geometrical kink shape')


sigma_b = array([0.001*sigma_p, 0.01*sigma_p, 0.1*sigma_p, 0.2*sigma_p, 0.25*sigma_p, \
                 0.3*sigma_p, 0.4*sigma_p, \
                 0.5*sigma_p, 0.6*sigma_p, 0.7*sigma_p, \
                 0.75*sigma_p, 0.8*sigma_p,  0.9*sigma_p, 0.99*sigma_p, 0.999*sigma_p])

(u_0, u_max, H_kp, Un, zdiff_kp) = kp_energy(xfine, yfine, yderfine, sigma_b)


plt.figure(4)
plt.plot(xfine,yderfine,'--',u_0,sigma_b,'o',u_max,sigma_b,'+',array([xfine[sigma_p_index]]),array([sigma_p]),'x')
plt.title('Derivative of potential')
plt.ylabel('dU/du')
plt.xlabel('u')


print     "Sigma*b (Pa)     u_0 (m)        u_max (m)   H_kp (kJ/mol) Un"
for i in xrange(len(sigma_b)):
    print "% .3g    % .3g      % .3g    % .3g    % .3g" % ((sigma_b[i]/5.0E-10), u_0[i], u_max[i], (H_kp[i]*Na/1000), (Un[i]*Na/1000))

plt.figure(6)
plt.plot((sigma_b/5.0E-10),(Un*Na/1000),'--')
plt.title('Critical stress / kink energy')
plt.xlabel('Tau* - Pa')
plt.ylabel('Un - kJ/mol')

plt.figure(7)
plt.plot((sigma_b/sigma_p),(Un/(kink_energy2*2.0)),'--',(sigma_b/sigma_p),(Un/(kink_energy2*2.0)),'o')
plt.title('Critical stress / kink energy')
plt.xlabel('Tau*/sigma_p')
plt.ylabel('Un/2Uk')

###plt.figure(5)
###for i in xrange(len(sigma_b)):
###    plt.plot(xfine[(u_0_index[i]+1):u_max_index[i]],(zdiff_kp[i,(u_0_index[i]+1):u_max_index[i]]),'-')
###   # plt.plot(xfine[(u_0_index[i]+1):u_max_index[i]],(5.0-zdiff_kp[i,(u_0_index[i]+1):u_max_index[i]]),'-')
###plt.xlabel('u')
###plt.ylabel('z-z0')
###plt.title('Kinki-pair shape')
plt.show()
