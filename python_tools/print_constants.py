#!/usr/bin/env python3

from my_work_env import *

z = 0.104325 
a = 1 / ( z + 1 )
print( 'a: %g, z: %g. Dcom: %g [g:%g], Dlum: %g [g:%g]'\
    %(\
    a,\
    z,\
    mycc.D_c(z),\
    mycc.D_c(z) / mycc.gadget_length_in_cm,\
    mycc.D_l(z),\
    mycc.D_l(z) / mycc.gadget_length_in_cm 
) )
