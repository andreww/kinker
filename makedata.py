#!/usr/bin/python

from numpy import *


xdata = array([0.01E-10, 0.5E-10, 1.0E-10, 1.5E-10, 2.0E-10, 2.5E-10, 3.0E-10, 3.5E-10, 4.0E-10, 4.5E-10, 4.99E-10])

G = 60.0E9

ydata = ((5.0E-10/(2.0*pi)) * 5.0E-10 * (0.01 * G) * (1 - cos((2 * pi * xdata) / 5.0E-10)))

ydata = (ydata - average(ydata)) + 0.5*G*2.5E-19 

for i in xrange(len(xdata)):
    print xdata[i],  ydata[i]

