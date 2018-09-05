#!/usr/bin/env python3

from scipy.special import beta as sbeta
from scipy.special import betainc as sbetainc
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c as scc

LightSpeed = scc * 100

def beta( x, a, b ):
    return sbeta(a,b) * sbetainc(a,b,x)

def n( a, q, c ):
    return c * np.power( q, 1-a ) / ( a-1 )

def e_div_c2( a, q, c ):
   #t = LightSpeed**2 * c / ( a-1 )
   t = c / ( a-1 )
   t = t * ( 0.5 * beta( 1/(1+q*q), (a-2)/2, (3-a)/2 ) + np.power( q, 1-a ) * ( np.sqrt( 1+q*q ) - 1 ) )
   return t

def p_div_c2( a, q, c ):
    #t = c * LightSpeed**2 / 6 * beta( 1/(1+q*q), (a-2)/2, (3-a)/3 )
    t = c / 6 * beta( 1/(1+q*q), (a-2)/2, (3-a)/2 )
    return t

def a2q2( a, A, B ):
    q2 = ( (A-3*B/(a-1)) + 1 ) ** 2 - 1
    return q2

def F( a, A, B ):
    q2 = a2q2( a, A, B )
    #print( a )
    r = ( a-1 ) / ( 6 * np.power(q2, (1-a)/2) ) * beta( 1/(1+q2), (a-2)/2, (3-a)/2 )
    r = B - r
    return r

fig = plt.figure()
#ax_AB = fig.add_subplot( 121 )
ax_F = fig.add_subplot( 211 )
ax_logF = fig.add_subplot( 212 )

aa = np.linspace( 2.1, 2.9, 100 )
aa2 = np.array( [2.2, 2.4, 2.7] )
qa2 = np.array( [10, 100, 10000] )
c = 1

for a in aa2:
    for q in qa2:

        nn = n( a, q, c )
        ee = e_div_c2( a, q, c )
        pp = p_div_c2( a, q, c )

        A = ee / nn
        B = pp / nn

        r = F( aa, A, B )

        qq = '%e'%(q)
        t = qq.split( 'e' )
        q1 = float( t[0] )
        q2 = int( t[1] )

        AA = '%e'%(A)
        t = AA.split( 'e' )
        A1 = float( t[0] )
        A2 = int( t[1] )

        BB = '%e'%(B)
        t = BB.split( 'e' )
        B1 = float( t[0] )
        B2 = int( t[1] )

        ax_F.plot( aa, r, '.', label=r'$\alpha=%.1f, q=%.1f \times 10^{%i}, A=%.1f \times 10^{%i}, B=%.1f \times 10^{%i}$'
                %(a, q1, q2, A1, A2, B1, B2) )
        ax_F.plot( aa, r )

        ax_logF.plot( aa, r, '.', label=r'$\alpha=%.1f, q=%.1f \times 10^{%i}, A=%.1f \times 10^{%i}, B=%.1f \times 10^{%i}$'
                %(a, q1, q2, A1, A2, B1, B2) )
        ax_logF.plot( aa, r )

ax_F.set_ylabel( r'$F$' )
ax_F.set_xlabel( r'$\alpha$' )
ax_F.legend()
ax_F.grid()

ax_logF.set_ylabel( r'$F$' )
ax_logF.set_xlabel( r'$\alpha$' )
ax_logF.set_yscale( 'log' )
ax_logF.legend()
ax_logF.grid()

plt.show()



