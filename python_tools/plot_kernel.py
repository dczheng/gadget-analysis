#!/usr/bin/env python3
from my_work_env import *
from scipy.integrate import quad

def k1( q ): # standar kernel
    norm = 8.0 / np.pi
    if 0 <= q <= 0.5:
        return norm * ( 1 - 6*q*q + 6*q*q*q )
    elif q <= 1:
        return norm * ( 2.0 * np.power( 1-q, 3.0 ) )
    else:
        return 0

def k2( q ): # quinitic kernel
    u1 = 1 - q
    u2 = 2/3 - q
    u3 = 1/3 - q
    norm = 2187 / ( 40 * np.pi )
    if 0<= q < 1/3:
        return (u1**5 - 6*u2**5 + 15*u3**5) * norm
    elif q < 2/3:
        return (u1**5 - 6*u2**5) * norm
    elif q <= 1:
        return (u1**5) * norm
    else:
        return 0

def k3( q ): #wendland c4 kernel
    norm = 495 / ( 32 * np.pi )
    u = 1-q
    if q<=1:
        return u**6 * ( 1 + 6*q + 35/3*q*q ) * norm
    else:
        return 0
def k4( q ): #wendland c6 kernel
    norm = 1365 / ( 64 * np.pi )
    u = 1-q
    if q<=1:
        return u**8 * ( 1 + 8*q + 25*q**2 + 32*q**3 ) * norm
    else:
        return 0

fig = plt.figure()
ax1 = fig.add_subplot( 221 )

ax2 = fig.add_subplot( 222 )

ax3 = fig.add_subplot( 223 )

ax4 = fig.add_subplot( 224 )


N = 100
q = np.linspace( 0, 1, 100 )
y1 = np.zeros( N )
y2 = np.zeros( N )
y3 = np.zeros( N )
y4 = np.zeros( N )
for i in range( N ):
    y1[i] = k1( q[i] )
    y2[i] = k2( q[i] )
    y3[i] = k3( q[i] )
    y4[i] = k4( q[i] )
ax1.plot( q, y1 )
ax1.set_title( "standard kernel" )
ax2.plot( q, y2 )
ax2.set_title( "quintic kernel" )
ax3.plot( q, y3 )
ax3.set_title( "wendland c4 kernel" )
ax4.plot( q, y4 )
ax4.set_title( "wendland c6 kernel" )

f1 = lambda x: k1(x) * x**2
f2 = lambda x: k2(x) * x**2
f3 = lambda x: k3(x) * x**2
f4 = lambda x: k4(x) * x**2

print( quad( f1, 0, 1 )[0] * np.pi * 4 )
print( quad( f2, 0, 1 )[0] * np.pi * 4 )
print( quad( f3, 0, 1 )[0] * np.pi * 4 )
print( quad( f4, 0, 1 )[0] * np.pi * 4 )

plt.savefig( 'kernel.pdf' )
