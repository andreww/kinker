#!/usr/bin/python

from numpy import *

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


def pot_func(x, p,x_len):
    N = (len(p)/2)
    a = p[0:N]
    b = p[N:2*N]
    x0 = 0.0
    value = 0
    for i in range(len(a)):
        value=value + (a[i]**4)*sin((x/x_len)*pi/2.0) \
                               *((sin((x/x_len)*pi))**b[i]) \
                               *cos((x/x_len)*pi/2.0)
    value = value + x0
    return value

def pot_num_deriv(x,p,x_len,delta_x):
    # Could do this analytically I think.
    y_plus = pot_func((x+delta_x),p,x_len)
    y_mins = pot_func((x-delta_x),p,x_len)
    dyBydx = (y_plus-y_mins)/(2.0*delta_x)
    return dyBydx

def pot_error(p, x, y, x_len):
    error = pot_func(x,p,x_len) - y
    return error

def fit_pot(p0, x, y, x_len):
    from scipy import optimize
    opt_p, cov_x, infodict, mesg, s = optimize.leastsq(pot_error, p0, args=(x, y, x_len), full_output=1, maxfev=0)
    print mesg
    return opt_p
    


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    print "Testing the Peiels potential fitter..." 
    basename = raw_input("Basename for potetial input file (input file is basname.dat): ")
    filename = basename + '.dat'
    (x,y) = load_data(filename)

    p0 = [2.4E-10, 2.4E-10, 1, 2, 1, 2, 3, 4, 1, 2]
    opt_p = fit_pot(p0, x, y, max(x))
    

    print opt_p
    opt_y = pot_func(x,opt_p,max(x))
    print opt_y

    plt.figure(1)
    plt.plot(x,y,'o',x,opt_y,'-')
    plt.ylabel('U')
    plt.xlabel('u')
    plt.title('U(u)')
    plt.show()

    
