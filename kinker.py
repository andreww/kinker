#!/usr/bin/python

from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate
from peierls import kp_energy, num_integ, kinker
import peierls_pot as pot
from convert import *

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

xmax_input = raw_input("Provide a maximum value of x in m (empty string to use input file and 5.0E-10 for old results):")
if (xmax_input == ""):
    xmax_input = None
else:
    xmax_input = float(xmax_input)

Na = 6.022E23 # Avagadro's number

# Load data set, x is the displacment (u) and y is the energy (U)
#(x,y) = load_data(None) # For test data
(x,y) = load_data(filename)

(xfine, yfine, sigma_p, zdiff, u_0, u_max, H_kp, Un, \
 zdiff_kp, yderfine, sigma_b, sigma_p_index, kink_energy2) = kinker(x,y, maxx=xmax_input)

# Now do it again using a function...

func = pot.choose_func()

# Fit this data to a sin function
print "Fitting initial potential to sin function"
#p0 = [2.4E-10, 2.4E-10, 1, 2, 1, 2, 3, 4, 1, 2]
p0 = [2.4, 0.001, 0.001, 2.0, 2.0, 2.0]
opt_p = pot.fit_pot(p0, x, y, max(x), func)

# Initial solution with this function
(xfine_init, yfine, sigma_p, zdiff_init, u_0_init, u_max, H_kp, Un_init, \
 zdiff_kp, yderfine_init, sigma_b_init, sigma_p_index, kink_energy2_init) = kinker(x,y,G=60E9,method=func,params=opt_p)


# Load experimental data
exptdata = raw_input("Expt data (input is basname.dat): ")
(T_n,expt_tau) = load_data(exptdata)
# Convert from MPa to Pa and subtract tau_crit (15 MPa, T_crit is 500 K, handled in error func.)
tau_n = expt_tau*1E6
t_crit = float(raw_input("Estimate t_crit (where tau ~10 MPa, used to use 650K): "))
tau_crit = 10.0E6

# Plot data 
plt.figure(1)
plt.plot(m_2_Ang(x),Jpm_2_eVpAng(y),'o',m_2_Ang(xfine),Jpm_2_eVpAng(func(xfine,opt_p, max(x))),'-')
plt.ylabel('E (eV/Ang)')
plt.xlabel('x (Ang)')
plt.title('E(x)')

plt.figure(2)
plt.plot(m_2_Ang(xfine[1:]),m_2_Ang(zdiff),'-')
plt.ylabel('z (Ang)')
plt.xlabel('x (Ang)')
plt.title('Geometrical kink shape')

plt.figure(3)
plt.plot(m_2_Ang(xfine),Pa_2_MPa(yderfine),'--',m_2_Ang(u_0_init),Pa_2_MPa(sigma_b_init),'o') #,u_max,sigma_b,'+') # ,array([xfine[sigma_p_index]]),array([sigma_p]),'x')
plt.title('Derivative of potential')
plt.ylabel('dE/dx (MPa)')
plt.xlabel('x (Ang)')

plt.figure(4)
plt.plot((sigma_b_init/sigma_p),(Un_init/(kink_energy2_init*2.0)),'--',(sigma_b_init/sigma_p),(Un_init/(kink_energy2_init*2.0)),'o')
plt.title('CRSS / kink energy (relative)')
plt.xlabel('sigma*/sigma_p')
plt.ylabel('Un/2Uk')

plt.figure(5)
plt.plot(((Un_init/(kink_energy2_init*2.0))*t_crit),Pa_2_MPa(sigma_b_init+tau_crit),'--',T_n,Pa_2_MPa(tau_n),'o',)
plt.legend(('model', 'expt data'))
plt.title('CRSS / kink energy (absolute)')
plt.ylabel('CRSS (MPa)')
plt.xlabel('T (K)')


plt.show()
