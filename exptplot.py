#!/usr/bin/python

from numpy import *
import matplotlib.pyplot as plt

FILE = open('MgO_experiments.dat', 'r')
T,tau = loadtxt(FILE, unpack=True)
FILE.close()

plt.figure(1)
plt.plot(T,tau,'o')
plt.xlabel('T')
plt.ylabel('tau')
plt.title('crss')

T_n = T/750.0
tau_n = (tau-10.0)/55.0
plt.figure(2)
plt.plot(tau_n,T_n,'o')
plt.xlabel('tau_n')
plt.ylabel('T_n')
plt.title('crss')
plt.show()
