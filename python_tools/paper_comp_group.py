#!/usr/bin/env python3

from my_work_env import *
import pandas as pd
from scipy.optimize import curve_fit

def fit_f( x, a, b ):
    return x*a+b


g_dir = sys.argv[1]
g_cre_dir = sys.argv[2]
snapindex = int( sys.argv[3] )
fn_out = sys.argv[4]
g =    pd.read_csv( g_dir + "/group_%03i.csv"%snapindex  )
g_cre = pd.read_csv( g_cre_dir + "/group_%03i.csv"%snapindex  )

L = 100000
Lh = L / 2.0
g_m, g_n = g.shape
g_cre_m, g_cre_n = g_cre.shape
print( 'g: %i, %i'%(g_m, g_n) )
print( 'g_cre: %i, %i'%(g_cre_m, g_cre_n) )
print( 'BoxSize: %g'%L )

fs = 5
rows = 1
cols = 2
fig, axs = plt.subplots( rows, cols, figsize=(fs*cols, fs*rows) )
ax_h = axs[0]
ax_f  = axs[1]
x = g['vr'].to_numpy()
y = g_cre['vr'].to_numpy()
x = x[ np.abs(x) < 0.8 ]
y = y[ np.abs(y) < 0.8 ]
x = np.abs(x+0.5)
y = np.abs(y+0.5)
vmin = min( [ x.min(), y.min() ] )
vmax = max( [ x.max(), x.max() ] )
print( "vmin: %g, vmax: %g"%( vmin, vmax ) )
bins = np.linspace( vmin, vmax, 15 )
print( bins )
ax_h.plot( x, '.', label='sim' )
ax_h.plot( y, '*', label='cre' )
#ax_h.hist( x, bins=bins, label='sim' )
#ax_h.hist( y, bins=bins, label='cre', histtype='step' )
ax_h.legend()
#ax_h.set_yscale( 'log' )

xb = np.histogram( x, bins=bins )[0]
yb = np.histogram( y, bins=bins )[0]
print( xb, yb )

#index = (xb>0) * (yb>0)
#xb = xb[index]
#yb = yb[index]
#xb = np.log10(xb)
#yb = np.log10(yb)
r = curve_fit( fit_f, xb, yb )
yy = fit_f( xb, r[0][0], r[0][1] )
print( r )
ax_f.plot( xb, yb, '*' )
ax_f.plot( xb, yy )
fig.savefig( fn_out )
exit()


test_keys = [ 'cre_i', 'i', 'mtot', 'mgas', 'mdm', 'mstar', 'ek', 'ep', 'vr', 'r200' ]
test_keys2 = [ 'mtot', 'ek', 'ep', 'vr' ]
for k in test_keys:
    print( "%6s "%(k), end='' )
for k in test_keys2:
    if k == 'vr':
        print( "%5.2s "%(k), end='' )
        continue
    print( "%11s "%(k), end='' )
print()

plot_group = 0
fs = 2
N = g_cre_m 
cols = 5
rows = 6

if plot_group:
    fig = plt.figure( figsize=( fs*cols, fs*rows  ) )
    fig_cre = plt.figure( figsize=( fs*cols, fs*rows  ) )
    dx = 1/cols
    dy = 1/rows
    axs = []
    axs_cre = []
    for i in range( rows ):
        for j in range( cols ):
            axs.append(fig.add_axes( [j*dx, (rows-1-i)*dy, dx, dy] ))
            axs_cre.append(fig_cre.add_axes( [j*dx, (rows-1-i)*dy, dx, dy] ))

idxs = []
idxs_cre = []
for i in range( N ):
    flag = 0
    for j in range( g_m ):
        r = 0
        for k in [ 'x', 'y', 'z' ]:
            t = g_cre.loc[i, k] - g.loc[j, k]
            if t > Lh:
                t -= L
            if t < -Lh:
                t += L
            r += t**2
        r = r**0.5
        if r < g_cre.loc[i, 'r200']:
            flag += 1
            idx = j
    if flag != 1:
        print( "flag: %i, i: %i"%(flag, i) )
        print( g_cre.loc[i] )
    else:
        print( '%6i %6i '%(i, idx), end='' )
        for k in test_keys[2:]:
            print( '%6.2f '%( \
                (g_cre.loc[i,k]-g.loc[idx,k])/g.loc[idx,k] * 100 ),\
            end='')
        for k in test_keys2:
            if k == 'vr':
                print( "%5.2f "%g_cre.loc[i,k], end='' )
                continue
            print( '%11.3e '%g_cre.loc[i,k], end='')
        print()

        idxs.append( idx )
        idxs_cre.append( i )

for i in range( len(idxs) ):
        if not(plot_group):
            continue
        d_cre = np.loadtxt( g_cre_dir +\
                '/Density/Density_%03i_%04i_y.dat'%(snapindex, idxs_cre[i]) )
        d = np.loadtxt( g_dir + \
                '/Density/Density_%03i_%04i_y.dat'%(snapindex, idxs[i]) )
        ax = axs[i]
        ax_cre = axs_cre[i]

        L_cre = np.abs(d_cre[0,1]) * 2 / 1000
        L = np.abs(d[0,1]) * 2 / 1000

        m_tot_cre= d_cre[0,8:14].sum()
        m_tot = d[0,8:14].sum()

        d_cre = d_cre[1:,:]
        d = d[1:,:]

        ax.imshow( d, norm=mplc.LogNorm(), cmap=cm.hot )
        ax_cre.imshow( d_cre, norm=mplc.LogNorm(), cmap=cm.hot )

        a, b = d.shape
        ax.text( 0.03*b, 0.03*a, r'$L:%.2f \, h^{-1}{\rm Mpc}$'%(L) )
        ax_cre.text( 0.03*b, 0.03*a, r'$L:%.2f \, h^{-1}{\rm Mpc}$'%(L_cre) )

        mx = real2tex( m_tot, b=10, n=2 )
        ax.text( 0.03*b, 0.9*a, r'$M: %s\, h^{-1} M_{\odot}$'%(mx) )
        mx = real2tex( m_tot_cre, b=10, n=2 )
        ax_cre.text( 0.03*b, 0.9*a, r'$M: %s\, h^{-1} M_{\odot}$'%(mx) )

for i in range(rows*cols):
    if not(plot_group):
        continue
    axs[i].set_yticks( [] )
    axs_cre[i].set_xticks( [] )
    axs_cre[i].set_yticks( [] )
    axs[i].invert_yaxis()
    axs_cre[i].invert_yaxis()

if plot_group:
    fig.savefig( "g_%03i.pdf"%snapindex )
    fig_cre.savefig( "g_cre_%03i.pdf"%snapindex )

rows = 3
cols = 3
fs = 5
fig, axs = plt.subplots( rows, cols, figsize=(fs*cols, fs*rows) )
ax_gas = axs[0,0]
ax_star = axs[0,1]
ax_vr = axs[0,2]
ax_dif = axs[1,0]
ax_war = axs[1,1]
ax_hot = axs[1,2]
ax_con = axs[2,0]

def myplot( x, y, ax, t ):
    r = curve_fit( fit_f, x, y )
    yy = fit_f( x, r[0][0], r[0][1] )
    print( r[0] )
    ax.plot( x, y, 'r.' )
    ax.plot( x, yy, 'b-', label=r'$k=%.3f, b=%.3f$'%(r[0][0],r[0][1] ) )
    ax.plot( x, x, 'y-' )
    ax.legend()
    ax.set_title( t )

def myplot2( x, y, ax, t ):
    x = x[np.abs(x)<1]
    y = y[np.abs(y)<1]
    xmin = min( [x.min(), y.min()] )
    xmax = max( [x.max(), y.max()] )
    bins = np.linspace( xmin, xmax, 10 )
    print( bins )
    #ax.hist( x, histtype='step', bins=bins )
    ax.hist( x, bins=bins )
    ax.hist( y, histtype='step', bins=bins )
    ax.set_title( t )

print( 'gas:' )
x = np.array( g.loc[idxs, 'mgas'] / g.loc[idxs, 'mtot'] )
y = np.array( g_cre.loc[idxs_cre, 'mgas'] / g_cre.loc[idxs_cre, 'mtot'] )
x = g['mgas'] / g['mtot']
y = g_cre['mgas'] / g_cre['mtot']
myplot2( x, y, ax_gas, 'gas' )

print( 'star:' )
x = np.array( g.loc[idxs, 'mstar'] / g.loc[idxs, 'mtot'] )
y = np.array( g_cre.loc[idxs_cre, 'mstar'] / g_cre.loc[idxs_cre, 'mtot'] )
x = g['mstar'] / g['mtot']
y = g_cre['mstar'] / g_cre['mtot']
myplot2( x, y, ax_star, 'star' )

print( 'vr:' )
x = np.array( g.loc[idxs, 'vr'] )
y = np.array( g_cre.loc[idxs_cre, 'vr'] )
#x = np.abs( x+0.5 )
#y = np.abs( y+0.5 )
myplot( x, y, ax_vr, 'vr' )
#myplot2( g['vr'], g_cre['vr'], ax_vr, 'vr' )

print( "diffuse:" )
x = np.array( g.loc[idxs, 'mdiffuse'] / g.loc[idxs, 'mgas'] )
y = np.array( g_cre.loc[idxs_cre, 'mdiffuse'] / g_cre.loc[idxs_cre, 'mgas'] )
myplot( x, y, ax_dif, 'diffuse' )

print( "warmhot:" )
x = np.array( g.loc[idxs, 'mwarmhot'] / g.loc[idxs, 'mgas'] )
y = np.array( g_cre.loc[idxs_cre, 'mwarmhot'] / g_cre.loc[idxs_cre, 'mgas'] )
myplot( x, y, ax_war, 'warmhot' )

print( "hot:" )
x = np.array( g.loc[idxs, 'mhot'] / g.loc[idxs, 'mgas'] )
y = np.array( g_cre.loc[idxs_cre, 'mhot'] / g_cre.loc[idxs_cre, 'mgas'] )
myplot( x, y, ax_hot, 'hot' )

print( "condensed:" )
x = np.array( g.loc[idxs, 'mcondensed'] / g.loc[idxs, 'mgas'] )
y = np.array( g_cre.loc[idxs_cre, 'mcondensed'] / g_cre.loc[idxs_cre, 'mgas'] )
myplot( x, y, ax_con, 'condensed' )

fig.savefig( fn_out )

