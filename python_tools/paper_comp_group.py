#!/usr/bin/env python3

from my_work_env import *
import pandas as pd
from scipy.optimize import curve_fit

def fit_f( x, a, b ):
    return x*a+b


g = pd.read_csv( sys.argv[1] )
g_cre = pd.read_csv( sys.argv[2] )
fn_out = sys.argv[3]

L = 100000
Lh = L / 2.0
g_m, g_n = g.shape
g_cre_m, g_cre_n = g_cre.shape
print( 'g: %i, %i'%(g_m, g_n) )
print( 'g_cre: %i, %i'%(g_cre_m, g_cre_n) )
print( 'BoxSize: %g'%L )

test_keys = [ 'i', 'cre_i', 'r200', 'mtot', 'mgas', 'mdm', 'mstar', 'ek', 'ep', 'vr' ]
test_keys2 = [ 'mtot', 'ek', 'ep', 'vr' ]
for k in test_keys:
    print( "%6s "%(k), end='' )
for k in test_keys2:
    if k == 'vr':
        print( "%5.2s "%(k), end='' )
        continue
    print( "%11s "%(k), end='' )
print()

x_k = 'vr'
y_k = 'vr'
x_log = 0
y_log = 0
x = []
y = []
N = 25
for i in range( 25 ):
    flag = 0
    for j in range( g_cre_m ):
        r = 0
        for k in [ 'x', 'y', 'z' ]:
            t = g.loc[i, k] - g_cre.loc[j, k]
            if t > Lh:
                t -= L
            if t < -Lh:
                t += L
            r += t**2
        r = r**0.5
        if r < g.loc[i, 'r200']:
            flag += 1
            idx = j
    if flag != 1:
        print( "flag: %i, i: %i"%(flag, i) )
        print( g.loc[i] )
    else:
        print( '%6i %6i '%(i, idx), end='' )
        for k in test_keys[2:]:
            print( '%6.2f '%( \
                (g_cre.loc[idx,k]-g.loc[i,k])/g.loc[i,k] * 100 ),\
            end='')
        for k in test_keys2:
            if k == 'vr':
                print( "%5.2f "%g.loc[i,k], end='' )
                continue
            print( '%11.3e '%g.loc[i,k], end='')
        print()
        if i == 7:
            continue
        if i == 4:
            continue

        #x.append( g.loc[i,'mstar'] / g.loc[i, 'mtot'] )
        #y.append( g_cre.loc[idx,'mstar'] / g_cre.loc[idx, 'mtot']  )
        x.append( g.loc[i,x_k] )
        y.append( g_cre.loc[idx,y_k] )
        #x.append( (g.loc[i,x_k]+g_cre.loc[idx,x_k])/2 )
        #x.append( (g_cre.loc[idx,x_k]-g.loc[i,x_k])/g.loc[i,x_k] * 100 )
        #x.append( g_cre.loc[idx,x_k] )
        #y.append( (g_cre.loc[idx,y_k]-g.loc[i,y_k])/g.loc[i,y_k] * 100 )
        #y.append( g_cre.loc[idx,y_k] )
x = np.array( x )
y = np.array( y )
r = curve_fit( fit_f, x, y )
print( r[0][0] )
yy = fit_f( x, r[0][0], r[0][1]  )
if x_k == 'ep':
    x = np.abs(x)
plt.plot( x, y, '*' )
plt.plot( x, yy, label=r'$k=%.2f$'%(r[0][0]) )
plt.xlabel( x_k )
plt.ylabel( y_k + '-diff' )
plt.legend()
if x_log:
    plt.xscale( 'log' )
if y_log:
    plt.yscale( 'log' )
plt.savefig( fn_out )
