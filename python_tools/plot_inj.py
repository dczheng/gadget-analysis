#!/usr/bin/env python3

import matplotlib as mpl
mpl.use( 'agg' )
import matplotlib.pyplot as plt
import numpy as np
import sys

data = np.loadtxt( './inj_test.dat' )

plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )


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


fig = plt.figure( figsize=(5, 5) )
ax_i = 1

index = 0
while( index < L ):
    '''
    if index != 10 and index != 43:
        index += 1
        continue
    '''
    if index != 10:
        index += 1
        continue

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

    #ax = fig.add_subplot( 120 + ax_i )
    #ax_i += 1
    ax = fig.add_subplot( 110 + ax_i )

    for j in range(3):
        c = t[j*4+0]
        a = t[j*4+1]
        qmin = t[j*4+2]
        qmax = t[j*4+3]
        f = calc_f( c, a, qmin, qmax, q )
        qmax_str = "%.1e"%(qmax)
        qmax_str = qmax_str.split( 'e' )
        #print( qmax_str[0], qmax_str[1][2:] )
        if j == 0:
            #ax.plot( q, f, line_style[j], label=r'$f_{\rm ori}\;\tilde{C}:%.1f,\;\alpha:%.1f,\;p_{\rm min}: %.1f,\; p_{\rm max}: %s \times 10^{%s}$'%(c, a, qmin, qmax_str[0], qmax_str[1][2:]) )
            ax.plot( q, f, line_style[j], label=r'$f_{\rm ori}$' )
        if j == 1:
            #ax.plot( q, f, line_style[j], label=r'$f_{inj}\;\tilde{C}:%.1f,\;\alpha:%.1f,\;p_{\rm min}: %.1f,\; p_{\rm max}: %s \times 10^{%s}$'%(c, a, qmin, qmax_str[0], qmax_str[1][2:]) )
            ax.plot( q, f, line_style[j], label=r'$f_{\rm inj}$' )
        if j == 2:
            #ax.plot( q, f, line_style[j], label=r'$f_{mod}\;\tilde{C}:%.1f,\;\alpha:%.1f,\;p_{\rm min}: %.1f,\; p_{\rm max}: %s \times 10^{%s}$'%(c, a, qmin, qmax_str[0], qmax_str[1][2:]) )
            ax.plot( q, f, line_style[j], label=r'$f_{\rm mod}$' )


    f1f2 = np.zeros( len(q) )

    for j in range(2):
        c = t[j*4+0]
        a = t[j*4+1]
        qmin = t[j*4+2]
        qmax = t[j*4+3]
        f = calc_f( c, a, qmin, qmax, q )
        f1f2 = f1f2 + f

    ax.plot( q, f1f2, line_style[3], label=r'$f_{\rm ori}+f_{\rm inj}$' )



    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    ax.set_xlabel( 'p', fontsize=20 )
    ax.set_ylabel( 'f', fontsize=20 )
    ax.yaxis.set_major_locator(plt.LogLocator(numticks=6))
    ax.minorticks_off()
    ax.tick_params( axis='both', labelsize=20, direction='in', pad=5 )
    ax.grid()
    ax.legend(prop={'size':15})
    index += 1


plt.tight_layout()
plt.savefig( 'inj.pdf' )

