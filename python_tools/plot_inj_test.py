#!/usr/bin/env python3

import matplotlib as mpl
mpl.use( 'agg' )
import matplotlib.pyplot as plt
import numpy as np
import sys

data = np.loadtxt( './inj_test.dat' )


def calc_f( c, a, qmin, qmax, pa ):
    r = []
    for p in pa:
        if p<qmin or p>qmax:
            r.append( 0 )
        else:
            r.append(c * p**(-a))

    #print( r )

    return np.array(r)

L, t = data.shape


index = 0
while( index < L ):
    t = data[index]

    qmin = 1e20
    qmax = -1

    for j in range(3):
        if t[j*4+2] < qmin:
            qmin = t[j*4+2]
        if t[j*4+3] > qmax:
            qmax = t[j*4+3]
    qmin = qmin - 1
    qmax = qmax + 1
    if ( qmin < 0 ):
        qmin = 1e-3
    print( "qmin: %g, qmax: %g"%(qmin, qmax) )

    q =  np.linspace( np.log10(qmin), np.log10(qmax), 30 )
    q = np.power( 10, q )

    line_style = [ '-', '-', '*', '-.' ]

    for j in range(3):
        c = t[j*4+0]
        a = t[j*4+1]
        qmin = t[j*4+2]
        qmax = t[j*4+3]
        f = calc_f( c, a, qmin, qmax, q )
        qmax_str = "%.1e"%(qmax)
        qmax_str = qmax_str.split( 'e' )
        #print( qmax_str[0], qmax_str[1][2:] )
        plt.plot( q, f, line_style[j], label=r'$f%i\;c=%.1f,\;\alpha=%.1f,\;q=(%.1f, %s \times 10^{%s})$'%(j+1, c, a, qmin, qmax_str[0], qmax_str[1][2:]) )

    f1f2 = np.zeros( len(q) )

    for j in range(2):
        c = t[j*4+0]
        a = t[j*4+1]
        qmin = t[j*4+2]
        qmax = t[j*4+3]
        f = calc_f( c, a, qmin, qmax, q )
        f1f2 = f1f2 + f
    plt.plot( q, f1f2, line_style[3], label=r'$f1+f2$' )



    plt.xscale( 'log' )
    plt.yscale( 'log' )
    plt.xlabel( 'q' )
    plt.ylabel( 'f' )
    plt.legend()

    plt.savefig( 'inj_test_eps/%i.eps'%(index) )
    #plt.savefig( 'inj_test/%i.png'%(index) )
    plt.close()

    index += 1
