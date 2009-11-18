#!/usr/bin/python

from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate, integrate

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

def num_integ(inx,iny):
    """Given values of a function (iny) evaluated at points (inx), 
       calculate integral of function at points using trapezium rule"""
    inx = atleast_1d(inx)
    out = zeros(inx.shape, dtype=inx.dtype)
    for n in xrange(len(out)):
        out[n] = integrate.trapz(iny[0:n],inx[0:n])
    return out


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


# Lower limits (u_0) come from eq.4:

sigma_b = array([0.001*sigma_p, 0.01*sigma_p, 0.1*sigma_p, 0.25*sigma_p, \
                 0.5*sigma_p, 0.75*sigma_p, 0.9*sigma_p, 0.99*sigma_p, 0.999*sigma_p])
u_0 = zeros(sigma_b.shape, dtype=sigma_b.dtype)
u_0_index = zeros(sigma_b.shape, dtype=integer) # What index do we want to start from?
for i in xrange(len(sigma_b)):
    for n in xrange(len(xfine)):
        if (yderfine[n] > sigma_b[i]):
            u_0[i] = xfine[n]
            u_0_index[i] = n
            break

# Upper limits (u_max) come from eq.9: this is a brute force search.

u_max = zeros(sigma_b.shape, dtype=sigma_b.dtype)
u_max_index = zeros(sigma_b.shape, dtype=integer)
U_u_max = zeros(sigma_b.shape, dtype=sigma_b.dtype)
inner_func = zeros((len(sigma_b),len(xfine)))

for i in xrange(len(sigma_b)):
    # Build an array inner_func = U(u) - U(u0) - (u - u_0)b_sigma
    # First index over values of b_sigma and second over values of u.
    # This is the inner bit of equation 8 and the LHS of equation 9.
    inner_func[i,:] = (yfine - yfine[u_0_index[i]] - ((xfine-u_0[i])*sigma_b[i]))

    for n in xrange(len(xfine)):
        if (n < (u_0_index[i]+5)):
            continue # u_max must be after u_0 and sometimes these seems to be a numerical
                     # issue so make it at least 5 spaces after u_0
        if (inner_func[i,n] < 0.0):
            u_max[i] = xfine[n]
            u_max_index[i] = n
            U_u_max[i] = yfine[n]
            break


plt.figure(4)
plt.plot(xfine,yderfine,'--',u_0,sigma_b,'o',array([xfine[sigma_p_index]]),array([sigma_p]),'x')
plt.title('Derivative of potential')
plt.ylabel('dU/du')
plt.xlabel('u')

# Calculate kink pair formation v's stress and report

H_kp = zeros(sigma_b.shape, dtype=sigma_b.dtype)
for i in xrange(len(sigma_b)):
    H_kp[i] = 2.0 * sqrt(2.0*Sd) * \
                (num_integ(    \
                  xfine[(u_0_index[i]):u_max_index[i]], \
                  sqrt(inner_func[i,(u_0_index[i]):u_max_index[i]]))[-1])
print     "Sigma*b (Pa)     u_0 (m)        u_max (m)   H_kp (kJ/mol)"
for i in xrange(len(sigma_b)):
    print "% .3g    % .3g      % .3g    % .3g" % ((sigma_b[i]/5.0E-10), u_0[i], u_max[i], (H_kp[i]*Na/1000))

# Calculate shape of kink pair - this looks odd? 

zdiff_kp = zeros((len(sigma_b),len(xfine)))
for i in xrange(len(sigma_b)):
    zdiff_kp[i,(u_0_index[i]+1):u_max_index[i]] = (sqrt(Sd/2.0) * \
                 (num_integ(xfine[(u_0_index[i]+1):u_max_index[i]], \
                 (1.0/(sqrt(yfine[(u_0_index[i]+1):u_max_index[i]]-yfine[u_0_index[i]]))))))


plt.figure(5)
for i in xrange(len(sigma_b)):
    plt.plot(xfine[(u_0_index[i]+1):u_max_index[i]],(zdiff_kp[i,(u_0_index[i]+1):u_max_index[i]]),'-')
   # plt.plot(xfine[(u_0_index[i]+1):u_max_index[i]],(5.0-zdiff_kp[i,(u_0_index[i]+1):u_max_index[i]]),'-')
plt.xlabel('u')
plt.ylabel('z-z0')
plt.title('Kinki-pair shape')
plt.show()
