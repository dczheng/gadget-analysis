#!/usr/bin/env python3

import matplotlib as mpl
mpl.use( 'agg' )
import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt( './cooling_test.dat' )
q = data[0, 6:]
f = data[1, 6:]

data = data[2:, :]

m, n = data.shape

mm = m // 2

for i in range( mm ):

    df_i = i * 2
    ff_i = i * 2 + 1

    rho   = data[df_i, 0]
    B     = data[df_i, 1]
    qmin  = data[df_i, 2]
    fqmin = data[df_i, 3]
    qmax  = data[df_i, 4]
    fqmax = data[df_i, 5]
    df    = data[df_i, 6:]
    ff    = data[ff_i, 6:]

    plt.plot( q, f, label=r'$f_t(p)$' )
    plt.plot( q, -df, '-.', label=r'$-df_t$' )
    plt.plot( q, f+df, '--', label=r'$f_t+df_t$' )
    plt.plot( q, ff, '.', label=r'$f_{t+\Delta t}$' )
    plt.plot( [qmin], [fqmin], '*', markersize=10, label=r'$ q_{min}$' )
    plt.plot( [qmax], [fqmax], '*', markersize=10, label=r'$ q_{max}$' )

    rho = "%.2e"%rho
    rho = rho.split( 'e' )

    if rho[1][0] == '-':
        t = '-'
    else:
        t = ''

    a = 1
    while( rho[1][a] == '0' ):
        a += 1
    t = t + rho[1][a:]
    rho[1] = t

    B = "%.2e"%B
    B = B.split( 'e' )
    if ( B[1][0] == '-' ):
        t = '-'
    else:
        t = ''

    a = 1
    while( B[1][a] == '0' ):
        a += 1
    t = t + B[1][a:]
    B[1] = t

    plt.text( 10, 1,
             r'$\frac{\rho}{\rho_{Bar}}: %s \times 10^{%s}\quad B: %s \times 10^{%s}$'%( rho[0],rho[1], B[0], B[1] ) )

    plt.xscale( 'log' )
    plt.yscale( 'log' )
    plt.xlabel( 'q' )

    plt.legend()
    plt.grid()

    #plt.savefig( 'cooling_test_%i.png'%i )
    plt.savefig( 'cooling_test_%i.eps'%i )
    plt.close()
