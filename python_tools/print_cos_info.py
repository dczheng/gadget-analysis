#!/usr/bin/env python3

from my_work_env import *

z = 0.01
while z<1:
    print( 
    "%.4f, %.4f,\
    %.2f, %.2f, %.2f,\
    %.2f, %.2f"%(z, 1/(1+z),\
    mycc.D_c(z)/mycc.Mpc,\
    mycc.D_a(z)/mycc.Mpc,\
    mycc.D_l(z)/mycc.Mpc,\
    mycc.D_a(z)/mycc.D_c(z), mycc.D_l(z)/mycc.D_c(z) ) \
    )
    z *= 1.1

print( "h/k_b:", mycc.h/mycc.k_b )
