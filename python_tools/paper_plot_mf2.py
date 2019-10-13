#!/usr/bin/env python3

from my_work_env import *
from scipy.optimize import curve_fit

def fit_f( x, a, b ):
    return x * a + b

d = np.loadtxt( sys.argv[1] )
d_cre = np.loadtxt( sys.argv[2] )
fn_out = sys.argv[3]
x = d[:,-1]
y = d_cre[:,-1]

i = 10

x = x[i:]
y = y[i:]
r = curve_fit( fit_f, x, y )
print( r[0] )
yy = fit_f( x, r[0][0], r[0][1] )

plt.plot( x, y, '.' )
plt.plot( x, yy, label=r'$a:%.2f, b:%.2f$'%( r[0][0], r[0][1] ) )
plt.xscale( 'log' )
plt.yscale( 'log' )
plt.legend()
plt.savefig( fn_out )

