#!/usr/bin/env python3

import matplotlib as mpl
#mpl.use( 'agg' )
from scipy.special import beta as sbeta
from scipy.special import betainc as sbetainc
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c as scc
from scipy.integrate import quad

LightSpeed = scc * 100

def beta( x, a, b ):

    if x<1e-30:
        return 0

    #print( x )

    if b>0:
        r = sbeta(a,b) * sbetainc(a,b,x)
        return r
    else:
        f = lambda x: x**(a-1) * (1-x)**(b-1)
        r = quad( f, 1e-20*x, x )
        return r[0]

def n( a, q, c ):
    return c * np.power( q, 1-a ) / ( a-1 )

def e_div_c2( a, q, c ):
   #t = LightSpeed**2 * c / ( a-1 )
    t = [ c / ( aa-1 ) * ( 0.5 * beta( 1/(1+q*q), (aa-2)/2, (3-aa)/2 ) + np.power( q, 1-aa ) * ( np.sqrt( 1+q*q ) - 1 ) ) for aa in a ]
    t = np.array( t )
    #t = c / (a-1) * ( 0.5 * b + np.power( q, 1-a ) * ( np.sqrt( 1+q*q ) - 1 ) )
    return t

def p_div_c2( a, q, c ):
    #t = c * LightSpeed**2 / 6 * beta( 1/(1+q*q), (a-2)/2, (3-a)/3 )
    t = [ c / 6.0 * beta( 1/(1+q*q), (aa-2)/2, (3-aa)/2 ) for aa in a ]
    t = np.array( t )
    return t

def a2q2( a, A, B ):
    q2 = ( (A-3*B/(a-1)) + 1 ) ** 2 - 1
    return q2

def F( a, A, B ):
    #print( a )

    r = [ (aa-1) / ( 6*np.power( a2q2(aa,A,B), (1-aa)/2 ) ) * beta( 1/(1+a2q2(aa,A,B)), (aa-2)/2, (3-aa)/2 ) for aa in a ]
    #r = ( a-1 ) / ( 6 * np.power(a2q2(a,A,B), (1-a)/2) ) * b
    r = B - r
    return r

def plot_test2():

    a = np.linspace( 2.1, 3.2, 10000 )

    c = 1
    aa = np.array([2.2])
    qq = 10

    nn = n( aa, qq, c )
    ee = e_div_c2( aa, qq, c )
    pp = p_div_c2( aa, qq, c )

    A = ee / nn
    B = pp / nn

    aa1 = np.array( [2.3] )
    qq1 = 10
    nn1 = n( aa1, qq1, c )
    ee1 = e_div_c2( aa1, qq1, c )
    pp1 = p_div_c2( aa1, qq1, c )

    A1 = ee1 / nn1
    B1 = pp1 / nn1

    A2 = ( ee+ee1 ) / (nn+nn1)
    B2 = ( pp+pp1 ) / ( nn+nn1 )

    print( "A: %f, B: %f"%(A,B) )
    print( "A1: %f, B1: %f"%(A1,B1) )
    print( "A2: %f, B2: %f"%(A2,B2) )

    r = F( a, A, B )
    r1 = F( a, A1, B1 )
    r2 = F( a, A2, B2 )

    fig = plt.figure()
    ax = fig.add_subplot( 111 )
    ax.plot( a, r, label='0' )
    ax.plot( a, r1, label='1' )
    ax.plot( a, r2, label='2' )
    plt.legend()
    ax.set_yscale( 'log' )
    plt.show()



def plot_test():
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

def value_test():
    c = 1
    a = 2.7
    q = 1000

    nn = n( a, q, c )
    ee = e_div_c2( a, q, c )
    pp = p_div_c2( a, q, c )
    A = ee / nn
    B = pp / nn

    a = 2.698
    q2 = a2q2( a, A, B )
    qq = np.sqrt( q2 )
    FF = F( a, A, B )
    bb  = beta( 1/(1+qq*qq), (a-2)/2, (3-a)/2 )
    x = (a-1) / ( 6 * np.power( qq, 1-a ) )

    print( " %.2f\n %.10e\n %.10e\n %.10e\n %.10e\n %.10e\n %.10e\n %.10e\n %.10e\n %.10e\n %.10e\n"%( a, FF, nn, ee, pp, A, B, x, bb, q2, qq ) )


#value_test()
#plot_test()
plot_test2()

