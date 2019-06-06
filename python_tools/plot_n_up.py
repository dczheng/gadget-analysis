#!/usr/bin/env python3


from my_work_env import *
import tools_and_constants as tc

def f(g, p, g_max, g_min):
    if g < g_max:
        if g > g_min:
            return np.power( g, -p ) * np.power( (1-g/g_max), p-2 )
        else:
            return np.power( g+g_min, -p )
    else:
        return 0

B = 3.24
z = 0.1
b1 =  3 * tc.sigma_t / ( 3 * tc.m_ec ) * 1e-12 / ( 8 * np.pi ) * B**2 * (1+z)**4
t = 10
print( 'b1:', b1 )
print( '1/b1: %f Myr'%( 1/b1 / tc.Myr) )
print( tc.Myr )

g_max = 1 / ( tc.Myr * b1 * t )
print( g_max )

g = np.power(10, np.linspace(3, np.log10(g_max), 1000) )
#g = np.power(10, np.linspace(2, 10, 1000) )
pa = [2.5, 3]
g_min = -1

for p in pa:
    N = [ f( i, p, g_max, g_min) for i in g ]
    print( p )
    #print( N )
    plt.loglog( g, N, label=r'$\alpha: %.1f$'%p )
    #plt.plot( g, N, label=r'$\alpha: %.1f$'%p )

plt.legend()
plt.minorticks_off()
plt.tick_params( axis='both', direction='in', labelsize=15, length=0 )
plt.xlabel( r'$\gamma$', fontsize=20 )
plt.ylabel( r'$N(\gamma, t)/N_0$', fontsize=20 )
plt.tight_layout()

plt.savefig( output_dir + 'n_up.pdf' )
