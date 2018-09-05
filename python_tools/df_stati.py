#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import random

def f( x, a, b ):
    return x**( a*x + b )

f1 = lambda x: f(x, -0.4, -1)

#print( quad( f1, 0, 10 ) )
N = 1000
x = np.linspace( 0.001, 10, N )
y1 = np.zeros( N )
y2 = np.zeros( N )
y3 = np.zeros( N )

for i in range( 10 ):
    a = -random.random()
    b = -random.random() * 2
    while a == 0:
        a = -random.random()
    print( a, b )
    f1 = lambda x: f(x, a, b)
    y1 = y1 + f1( x )

for i in range( 1000 ):
    a = -random.random()
    b = -random.random() * 2
    while a == 0:
        a = -random.random()
    print( a, b )
    f2 = lambda x: f(x, a, b)
    y2 = y2 + f2( x )

for i in range( 10000 ):
    a = -random.random()
    b = -random.random() * 2
    while a == 0:
        a = -random.random()
    print( a, b )
    f3 = lambda x: f(x, a, b)
    y3 = y3 + f3( x )



plt.plot( x, y1, '.', label='y1' )
plt.plot( x, y2, '.', label='y2' )
plt.plot( x, y3, '.', label='y3' )

plt.xscale( 'log' )
plt.yscale( 'log' )
plt.legend()
plt.show()
