#!/usr/bin/env python3

from my_work_env import *

fn2 = sys.argv[1]
fn2_cre = sys.argv[2]
fn0 = sys.argv[3]
fn0_cre = sys.argv[4]
fn_out = sys.argv[5]

fns = [ fn2, fn2_cre, fn0, fn0_cre ]
ds = []

for i in range( len(fns) ):
    print( 'load %s ...'%fns[i] )
    ds.append( np.loadtxt( fns[i] ) )

m, n = ds[0].shape
R = ds[0][1:,0] / 1000
N = ( n-1 ) // 3

if ( n-1 ) % 3 != 0:
    print( "error" )
    exit()

labels = []
for i in range( N ):
    t = real2tex( ds[0][0, 3*i+1], b=10 )
    if i == N-1:
        labels.append( r"$>%s M_{\odot}$"%t )
    else:
        t1 = real2tex( ds[0][0, 3*(i+1)+1], b=10 )
        labels.append( r"$%s \sim %sM_{\odot}$"%(t, t1) )
print( labels )
#labels_prefix = [ '', 'CRE\,', '', 'CRE\,' ]
cs = [ 'k', 'g', 'b', 'r' ]
ls = [ '-', '--', '-', '--' ]

fs = 5
r_diff= 0.3
margin = 0.15
fig = plt.figure( figsize=( fs*2/(1-margin), fs/(1-r_diff-margin) ) )
dx = (1-margin) / 2
dy = (1-margin-r_diff)

axs = []
for i in range(2):
    axs.append( fig.add_axes( [ margin/2+i*dx, r_diff+margin/2, dx, dy ] ) )
    axs.append( fig.add_axes( [ margin/2+i*dx, margin/2, dx, r_diff ] ) )

Nmin = 1000

Tunit = 1e7
vmin = 1e30
vmax = -vmin
for i in range(4):
    ax = axs[i//2*2]
    d = ds[i]
    for j in range(N):
        T = d[1:, j*3+1] / 1e7
        num = d[1:, j*3+3]
        if T.max() == 0:
            continue
        index = num > Nmin 
        if not( True in index ):
            continue

        x = R[index]
        y = T[index]
        ax.plot( x, y, cs[j%len(cs)] + ls[i], label=labels[j] )
        vmin = np.min( [T.min(), vmin] )
        vmax = np.max( [T.max(), vmax] )

print( vmin, vmax )
for i in range(4):
    ax = axs[i//2*2]
    ax.set_ylim( [vmin, vmax] )

evmin = 1e30
evmax = -evmin
for i in range(2):
    ax = axs[i*2+1]
    y = ds[i*2]
    y_cre = ds[i*2+1]
    for j in range(N):
        yy = y[1:, j*3+1] / Tunit
        yy_cre = y_cre[1:, j*3+1] / Tunit
        num = y[1:, j*3+3]
        num_cre = y_cre[1:, j*3+3] 
        #print ( yy.max(), yy_cre.max() )
        if yy.max() == 0 or yy_cre.max() == 0:
            continue

        index1 = yy > 0
        index2 = yy_cre > 0
        index3 = num > Nmin
        index4 = num_cre > Nmin 
        index = index1 * index2 * index3 * index4
        if not( True in index ):
            continue

        xx = R[index]
        yy = yy[index]
        yy_cre = yy_cre[index]
        yyerr = ( yy_cre - yy ) / yy
        #print( yyerr )
        #print( yyerr.min(), yyerr.max() )
        evmin = np.min( [yyerr.min(), evmin] )
        evmax = np.max( [yyerr.max(), evmax] )
        ax.plot( xx, yyerr, cs[j%len(cs)]+'-.', label=labels[j] )

print( evmin, evmax )
for i in range(2):
    ax = axs[i*2+1]
    ax.set_ylim( [evmin, evmax] )

for i in range(4):
    a = axs[i]
    a.legend( prop={'size':10}, framealpha=0.1  )
    a.set_xscale( 'log' )
    a.set_xlim( [R[0], R[-1]] )
    #a.grid()
    if i % 2 == 0:
        remove_tick_labels( a, 'x' )
    else:
        a.set_xlabel( r'$\rm Mpc$', fontsize=20 )

    '''
    if i == 0:
        a.set_title( r'$z=2$' )
    if i == 2:
        a.set_title( r'$z=0$' )
    '''

    if i > 1:
        remove_tick_labels( a, 'y' )
    else:
        a.set_ylabel( r'$\rm 10^7 \, K$', fontsize=20 )

    a.tick_params( axis='both', direction='in', labelsize=15 )


fig.savefig( fn_out )
