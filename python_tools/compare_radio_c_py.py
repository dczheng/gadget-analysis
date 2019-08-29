#!/usr/bin/env python3
    
from my_work_env import *
from radio import *

dd = np.loadtxt( sys.argv[1] )

m, n = dd.shape

for i in range(m):
    d = dd[i]
    nu = d[0]
    B = d[1]
    h = d[2]
    c = d[3]
    a = d[4]
    qmin = d[5]
    qmax = d[6]
    f_c = d[7]
    f_py = dp2dVdv( B, nu, c, a, qmin, qmax ) * ( 4/3 * np.pi * h**3 )
    print( "c: %10.3e, py: %10.3e, err: %6.2f%%"%( f_c, f_py, (f_c-f_py)/f_py*100 ) )


