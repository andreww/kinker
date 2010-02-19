#!/usr/bin/python

import numpy as np

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

def choose_func():
    i = 0
    funcs = (_pot_func,_simp_func) 
    for func in funcs:
        print str(i) + ": " + func.__doc__
        i = i + 1
    i = int(raw_input("Choose a function: "))
    return funcs[i]

def _pot_func(x, p,x_len):
    """Function of the form y = Sum{Ai*sin(x*pi/2)*sin^Bi(x/pi)*cos(x*pi/2)}"""
    N = (len(p)/2)
    a = p[0:N]
    b = p[N:2*N]
    x0 = 0.0
    value = 0
    for i in range(len(a)):
        value=value + (a[i]**4)*np.sin((x/x_len)*np.pi/2.0) \
                               *((np.sin((x/x_len)*np.pi))**(b[i]**2.0)) \
                               *np.cos((x/x_len)*np.pi/2.0)
    value = value + x0
    return value

def _simp_func(x, p,x_len):
    """Function of the form y = Sum{Ai*sin^Bi(x/pi)}"""
    N = (len(p)/2)
    a = p[0:N]
    b = p[N:2*N]
    x0 = 0.0
    value = 0
    for i in range(len(a)):
        value=value + (a[i]**4)*((np.sin((x/x_len)*np.pi))**(b[i]**2.0)) \
                               
    value = value + x0
    return value

def eval_func(x, p, x_len, func):
    return(func(x, p, x_len))

def pot_num_deriv(x,p,x_len,delta_x, func):
    # Could do this analytically I think.
    y_plus = func((x+delta_x),p,x_len)
    y_mins = func((x-delta_x),p,x_len)
    dyBydx = (y_plus-y_mins)/(2.0*delta_x)
    return dyBydx

def pot_error(p, x, y, x_len, func):
    error = func(x,p,x_len) - y
    return error

def fit_pot(p0, x, y, x_len, func):
    from scipy import optimize
    opt_p, cov_x, infodict, mesg, s = optimize.leastsq(pot_error, p0, args=(x, y, x_len, func), full_output=1, maxfev=0)
    print mesg
    return opt_p
    


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    print "Testing the Peiels potential fitter..." 

    func = choose_func()

    basename = raw_input("Basename for potetial input file (input file is basname.dat): ")
    filename = basename + '.dat'
    (x,y) = load_data(filename)

    p0 = [2.4E-10, 2.4E-10, 1, 2, 1, 2, 3, 4, 1, 2]
    opt_p = fit_pot(p0, x, y, max(x), func)
    

    print opt_p
    opt_y = eval_func(x,opt_p,max(x), func)
    print opt_y

    plt.figure(1)
    plt.plot(x,y,'o',x,opt_y,'-')
    plt.ylabel('U')
    plt.xlabel('u')
    plt.title('U(u)')
    plt.show()

    
