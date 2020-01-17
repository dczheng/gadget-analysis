#!/usr/bin/env python3

from my_work_env import *

z = float(sys.argv[1])
a = 1 / ( 1+z )

print( "z: %g"%z )
print( "a: %g"%a )
print( "Comoving distance %.2f Mpc"%( mycc.D_c(z)/mycc.Mpc) )
print( "Luminosity distance %.2f Mpc"%( mycc.D_l(z)/mycc.Mpc) )
print( "Angular diameter distance %.2f Mpc"%( mycc.D_a(z)/mycc.Mpc) )
print( "Myr: %g"%mycc.Myr )
print( "Mpc: %g"%mycc.Mpc )
print( mycc.gadget_gauss )
print( "1Myr * 1km/s = %g pc"%( 1e6 * mycc.yr * 1e5 / mycc.pc ) )
