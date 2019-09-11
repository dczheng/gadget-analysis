#!/usr/bin/python3
from os import listdir
import re

def check_completeness( f, s1, s2 ):
    lines = open( f ).readlines()
    l1 = []
    l2 = []
    for i in range(len(lines)):
        l = lines[i]
        #print( l )
        if s1 in l:
            l1.append( "line %i:  "%i + l )
        if s2 in l:
            l2.append( "line %i:  "%i + l )
    r = [ "%s: %i"%(s1, len(l1)), "%s: %i"%(s2, len(l2)) ]
    if len(l1) == len(l2):
        r.append( 'ok' )
        return r

    if l1 == [] and l2 != []:
        r.append(( ''.join(l2) ))
        return r

    #print( l1 )
    #print( l2 )
    flag = [0] * len(l1)
    for l in l2:
        i, j = re.search( '\(.*\)', l, re.DOTALL ).span()
        t = l[i+1:j-1].strip()
        #print( t )
        for j in range(len(l1)):
            if t in l1[j] and flag[j] == 0:
                flag[j] = 1
                break
    t = []
    for i in range(len(l1)):
        if flag[i] == 0:
            t.append( l1[i] )
    if r == []:
        r.append( 'ok' )
        return r
    else:
        r.append(( ''.join(t) ))
        return r

fs = listdir( '.' )
for f in fs:
    if f[-2:] != '.c':
        continue
    print( '-' * 50 )
    print( "%s ..."%f )
    r = check_completeness( f, 'fopen', 'fclose' )
    print( r[:-1] )
    print( 'check `fopen` and `fclose`: %s'%r[-1] )

    r = check_completeness( f, 'mymalloc', 'myfree' )
    print( r[:-1] )
    print( 'check `mymalloc` and `myfree`: %s'%r[-1] )
