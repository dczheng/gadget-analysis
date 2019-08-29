#!/usr/bin/env python3
    
from my_work_env import *

def f( c, a, q1, q2, q ):
    r = c * np.power( q, -a )
    r[ q<q1 ] = 0
    r[ q>q2 ] = 0
    return r

def compare_elecspec():

    N = 100
    fig, axs = plt.subplots( 1, 2, figsize=(2*5, 5) )
    for idx in range(fn_start_idx, fn_end_idx+1):
        fn = fn_pre + '%04i.dat'%idx
        print( "plot %s ..."%fn )
        d = np.loadtxt( fn )
        qmin = d[:,2][ d[:,2]>0 ].min() / 1.1
        qmax = d[:,3].max() * 1.1
        #print( qmin, qmax )
        q = np.logspace( np.log10(qmin), np.log10(qmax), N )
        F = np.zeros( [N] )
        
        hd = d[ 0, : ]
        d = d[1:,:]
        m, n = d.shape
        rho = 0
        for i in range( m ):
            dd = d[i,:]
            c = dd[0]
            a = dd[1]
            qmin = dd[2]
            qmax = dd[3]
            #print( dd )
            #r = f( dd[0], dd[1], dd[2], dd[3], q ) * dd[4]
            r = f( c, a, qmin, qmax, q )
            rho += dd[4]
            #print( r )
        
            F += r
    
        F = F / m
        axs[0].loglog( q, F, label="%i"%idx )
    axs[0].set_title( 'py' )
    
    d = np.loadtxt( fn_c )
    q = d[0,2:]
    for i in range( fn_start_idx+1, fn_end_idx+2 ):
        print( "plot gadget-tools's outputs [%i] ..."%(i-1) )
        F = d[i, 2:]
        axs[1].loglog( q, F, label="%i"%(i-1) )
    axs[1].set_title( 'c' )
    
    py_xlim = list(axs[0].get_xlim())
    c_xlim = list(axs[1].get_xlim())
    xlim = [ np.min( [py_xlim[0], c_xlim[0]] ), np.max( [py_xlim[1], c_xlim[1]] ) ]
    
    py_ylim = list(axs[0].get_ylim())
    c_ylim = list(axs[1].get_ylim())
    ylim = [ np.min( [py_ylim[0], c_ylim[0]] ), np.max( [py_ylim[1], c_ylim[1]] ) ]
    
    for i in range(2):
        axs[i].set_xlabel( 'q' )
        axs[i].set_ylabel( 'F' )
        axs[i].legend()
        axs[i].set_xlim( xlim )
        axs[i].set_ylim( ylim )
    
    
    plt.savefig( fn_out )

fn_pre = sys.argv[1]
fn_start_idx = int( sys.argv[2] )
fn_end_idx = int( sys.argv[3] )
fn_c = sys.argv[4]
fn_out = sys.argv[5]

compare_elecspec()
