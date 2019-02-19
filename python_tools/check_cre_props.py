#!/usr/bin/env python3

from elec_cooling import *

d = open( sys.argv[1] ).readlines()

def remove_comma( s ):
    if s[-1] == ',':
        return float( s[:-1] )
    else:
        return float( s )

tags = [ 'c:', 'a:', 'qmin:', 'qmax:', 'n:', 'e:' ]

for l in d:
    t = l[:-1].split()
    N = len( t )
    cre = []
    for i in range( N ):
        if t[i] in tags:
            cre.append( remove_comma( t[i+1] ) )
    if ( len(cre) != 6 ):
        continue

    if ( cre[0] == 0 ):
        print( '0' )
        continue

    n = cre_n2( cre[0], cre[1], cre[2], cre[3] )
    e1 = cre_e_gadget( cre[0], cre[1], cre[2] )
    e2 = cre_e_gadget_2( cre[0], cre[1], cre[2], cre[3] )

    #print( cre[4], n, cre[5], e )
    print( '%10s, qmin: %8.2f, qmax: %8.2e, n: %.4f%%, e1: %.4f%%(%.2g), e2: %.4f%%(%.2g)' \
            %( t[1], cre[2], cre[3], abs(n-cre[4])/n*100, \
                abs(e1-cre[5])/cre[5]*100, e1,  \
                abs(e2-cre[5])/cre[5]*100, e2 ) )
