#!/usr/bin/env python3

from my_work_env import *
import pandas as pd
from scipy.optimize import curve_fit

def fit_f( x, a, b ):
    return x*a+b


def comp():
   g_dir = sys.argv[1]
   g_cre_dir = sys.argv[2]
   snap_start = int( sys.argv[3] )
   snap_end = int( sys.argv[4] )
   fn_out = sys.argv[5]
   
   a = []
   b = []
   a_cre = []
   b_cre = []
   z = []
   r_gas = []
   r_gas_cre = []
   r_gas_var = []
   r_gas_cre_var = []
   rmeans = { 'mgas':[], 'mstar':[], 'mbaryon':[], 'T':[] }
   rmeans_cre = { 'mgas':[], 'mstar':[], 'mbaryon':[], 'T':[] }
   rvars = { 'mgas':[], 'mstar':[], 'mbaryon':[], 'T':[] }
   rvars_cre = { 'mgas':[], 'mstar':[], 'mbaryon':[], 'T':[] }
   
   keys = [ 'mgas', 'mstar', 'T' ]
   
   for snap_idx in range(snap_start, snap_end+1):
       fn     = "%s/group_%03i.csv"%(g_dir, snap_idx)
       fn_cre = "%s/group_%03i.csv"%(g_cre_dir, snap_idx)
   
       print( '-' * 20 )
       print( "idx: %i"%snap_idx )
   
       #g =    pd.read_csv( fn )
       #g_cre = pd.read_csv( fn_cre )
       gg =    pd.read_csv( fn )
       gg_cre = pd.read_csv( fn_cre )
   
       idx = gg['mtot'] <= 5e3
       g = gg[idx]
   
       idx = gg_cre['mtot'] <= 5e3
       g_cre = gg_cre[idx]
   
       print( g.shape, g_cre.shape )
   
       if g.shape[0] < 10 or g_cre.shape[0] <= 10:
           continue
   
       if g['z'].to_numpy()[0] != g_cre['z'].to_numpy()[0]:
           print( 'error' )
           exit()
   
       print( 'z:%g'%g['z'].to_numpy()[0] )
       print( "Ngroup: (%i, %i)"%(g.shape[0], g_cre.shape[0]) )
   
       for k in keys:
           if k == 'T':
               x = g[k] / 1e7
               y = g_cre[k] / 1e7
           else:
               x = g[k] / g['mtot']
               y = g_cre[k] / g_cre['mtot']
   
           rmeans[k].append( x.mean() )
           rvars[k].append( np.sqrt(x.var()) )
   
           rmeans_cre[k].append( y.mean() )
           rvars_cre[k].append( np.sqrt(y.var()) )
   
           print( '%s: %.3f[%.2e], %.3f[%.2e], (%.2f, %.2f) (%.2f, %.2f)'%( k, \
                   x.mean(), x.var(), y.mean(), y.var(), x.min(), x.max(), y.min(),\
                   y.max()
                   ) )
   
       x = (g['mstar'] + g['mgas']) / g['mtot']
       y = (g_cre['mstar'] + g_cre['mgas']) / g_cre['mtot']
       rmeans['mbaryon'].append( x.mean() )
       rvars['mbaryon'].append( np.sqrt(x.var()) )
       rmeans_cre['mbaryon'].append( y.mean() )
       rvars_cre['mbaryon'].append( np.sqrt(y.var()) )
       print( '%s: %.3f[%.2e], %.3f[%.2e]'%( 'mbaryon', \
                   x.mean(), x.var(), y.mean(), y.var() ) )
   
       z.append( g['z'].to_numpy()[0] )
   
   for k in rmeans.keys():
       rmeans[k] = np.array(rmeans[k])
       rmeans_cre[k] = np.array(rmeans_cre[k])
       rvars[k] = np.array(rvars[k])
       rvars_cre[k] = np.array(rvars_cre[k])
   
   
   fs = 5
   m = 1
   n = len(rmeans) 
   r_pad = 0.1
   r_err = 0.3
   
   dx = 1 / ( n + (n+1)*r_pad )
   dy = 1 / ( m + (m+1)*r_pad + r_err )
   
   fig = plt.figure( figsize=(fs/dx, fs/dy) )
   axs = [ fig.add_axes(\
       [ dx*(r_pad+j+j*r_pad), dy*(r_pad+r_err), dx, dy ] ) \
               for j in range(n) ]
   axs_err = [ fig.add_axes(\
       [ dx*(r_pad+j+j*r_pad), dy*r_pad, dx, dy*r_err ] ) \
               for j in range(n) ]
   
   keys.append( 'mbaryon' )
   cs = ['r', 'b']
   alpha = 0.5
   r_sigma = 3
   for k in range(len(keys)):
       ax = axs[k]
       ax_err = axs_err[k]
   
       ax.plot( z, rmeans[keys[k]], cs[0], label='sim' )
       ax.plot( z, rmeans_cre[keys[k]], cs[1], label='cre' )
   
       ax.fill_between( z,\
           rmeans[keys[k]] - r_sigma * rvars[keys[k]],\
           rmeans[keys[k]] + r_sigma * rvars[keys[k]],\
           alpha=alpha, facecolor=cs[0]\
       )
       ax.set_title( keys[k] )
   
       ax.fill_between( z,\
           rmeans_cre[keys[k]] - r_sigma * rvars_cre[keys[k]],\
           rmeans_cre[keys[k]] + r_sigma * rvars_cre[keys[k]],\
           alpha=alpha, facecolor=cs[1]\
       )
   
       ax_err.plot( z, rmeans_cre[keys[k]]-rmeans[keys[k]] )
   
       ax.legend()
   
   
   plt.savefig( fn_out )

def comp_gas():
    g_dir = sys.argv[1]
    g_cre_dir = sys.argv[2]
    snap_start = int( sys.argv[3] )
    snap_end = int( sys.argv[4] )
    fn_out = sys.argv[5]
    
    a = []
    b = []
    a_cre = []
    b_cre = []
    z = []
    r_gas = []
    r_gas_cre = []
    r_gas_var = []
    r_gas_cre_var = []
    rmeans = { 'mgas_r200':[], 'mcondensed':[], 'mhot':[], 'mdiffuse':[] }
    rmeans_cre = { 'mgas_r200':[], 'mcondensed':[], 'mhot':[], 'mdiffuse':[] }
    rvars = { 'mgas_r200':[], 'mcondensed':[], 'mhot':[], 'mdiffuse':[] }
    rvars_cre = { 'mgas_r200':[], 'mcondensed':[], 'mhot':[], 'mdiffuse':[] }
    
    keys = [ 'mgas_r200', 'mcondensed','mhot', 'mdiffuse' ]
    
    for snap_idx in range(snap_start, snap_end+1):
        fn     = "%s/group_%03i.csv"%(g_dir, snap_idx)
        fn_cre = "%s/group_%03i.csv"%(g_cre_dir, snap_idx)
        print( '-' * 20 )
        print( "idx: %i"%snap_idx )
    
        g =    pd.read_csv( fn )
        g_cre = pd.read_csv( fn_cre )
    
        if g['z'][0] != g_cre['z'][0]:
            print( 'error' )
            exit()
    
        print( 'z:%g'%g['z'][0] )
        print( "Ngroup: (%i, %i)"%(g.shape[0], g_cre.shape[0]) )
    
        for k in keys:
            '''
            if k == 'mgas_r200':
                x = g[k]
                y = g_cre[k]
            else:
                x = g[k] / g['mgas']
                y = g_cre[k] / g_cre['mgas']
            '''
            x = 1 - g[k] / g['mgas']
            y = 1 - g_cre[k] / g_cre['mgas']
            
    
            rmeans[k].append( x.mean() )
            rvars[k].append( np.sqrt(x.var()) )
    
            rmeans_cre[k].append( y.mean() )
            rvars_cre[k].append( np.sqrt(y.var()) )
    
            print( '%s: %.3f[%.2e], %.3f[%.2e]'%( k, \
                    x.mean(), np.sqrt(x.var()), y.mean(), np.sqrt(y.var()) ) )
    
        z.append( g['z'][0] )
    
    for k in rmeans.keys():
        rmeans[k] = np.array(rmeans[k])
        rmeans_cre[k] = np.array(rmeans_cre[k])
        rvars[k] = np.array(rvars[k])
        rvars_cre[k] = np.array(rvars_cre[k])
    
    
    fs = 5
    m = 1
    n = 4
    r_pad = 0.1
    r_err = 0.3
    
    dx = 1 / ( n + (n+1)*r_pad )
    dy = 1 / ( m + (m+1)*r_pad + r_err )
    
    fig = plt.figure( figsize=(fs/dx, fs/dy) )
    axs = [ fig.add_axes(\
        [ dx*(r_pad+j+j*r_pad), dy*(r_pad+r_err), dx, dy ] ) \
                for j in range(n) ]
    axs_err = [ fig.add_axes(\
        [ dx*(r_pad+j+j*r_pad), dy*r_pad, dx, dy*r_err ] ) \
                for j in range(n) ]
    
    cs = ['r', 'b']
    alpha = 0.5
    r_sigma = 1
    for k in range(len(keys)):
        ax = axs[k]
        ax_err = axs_err[k]
    
        ax.plot( z, rmeans[keys[k]], cs[0], label='sim' )
        ax.plot( z, rmeans_cre[keys[k]], cs[1], label='cre' )
    
        
        ax.set_title( keys[k].replace( '_', '-' ) )
    
        ax.fill_between( z,\
            rmeans[keys[k]] - r_sigma * rvars[keys[k]],\
            rmeans[keys[k]] + r_sigma * rvars[keys[k]],\
            alpha=alpha, facecolor=cs[0]\
        )

        ax.fill_between( z,\
            rmeans_cre[keys[k]] - r_sigma * rvars_cre[keys[k]],\
            rmeans_cre[keys[k]] + r_sigma * rvars_cre[keys[k]],\
            alpha=alpha, facecolor=cs[1]\
    )

        ax_err.plot( z, rmeans_cre[keys[k]]-rmeans[keys[k]] )

        ax.legend()

    plt.savefig( fn_out )

def find_match_id( gg1, gg2, m_min=0 ):

    idxs2 = []
    idxs1 = []
    idxs1_failed = []
    g1 = gg1[ gg1.mtot >= m_min ]
    g2 = gg2[ gg2.mtot >= m_min/1.5 ]

    g1_m, g1_n = g1.shape
    g2_m, g2_n = g2.shape
    L = 100000
    Lh = L / 2.0
    print( 'g1: %i, %i'%(g1_m, g1_n) )
    print( 'g2: %i, %i'%(g2_m, g2_n) )
    print( 'BoxSize: %g2'%L )
    for i in range( g1_m ):
        flag = 0
        for j in range( g2_m ):
            r = 0
            for k in [ 'x', 'y', 'z' ]:
                t = g1.loc[i, k] - g2.loc[j, k]
                if t > Lh:
                    t -= L
                if t < -Lh:
                    t += L
                r += t**2
            r = r**0.5
            if r < g1.loc[i, 'r200'] * 0.1:
                flag += 1
                idx = j
        if flag != 1:
            idxs1_failed.append(i)
            print( "can't match %i or over mach [%i]"%(i, flag) )
        else:
            idxs2.append( idx )
            idxs1.append( i )
    print( len(idxs1), len(idxs1_failed), len(idxs1_failed)/g2_m )
    return ( idxs1, idxs2, idxs1_failed )

def comp_single():
    g_dir = sys.argv[1]
    g_cre_dir = sys.argv[2]
    snapindex = int( sys.argv[3] )
    fn_out = sys.argv[4]
    g =    pd.read_csv( g_dir + "/group_%03i.csv"%snapindex  )
    g_cre = pd.read_csv( g_cre_dir + "/group_%03i.csv"%snapindex  )
    
    g_m, g_n = g.shape
    g_cre_m, g_cre_n = g_cre.shape
    
    idxs_cre, idxs, idxs_failed = find_match_id( g_cre, g, 1e4 )

    keys = [ 'mtot', 'mgas', 'mdm', 'mstar' ]
    key_fmts = {}
    for k in keys:
        key_fmts[k] = '%.2f,%.2f,%.2f,'
    #keys2 = ['mwarmhot', 'mhot', 'mcondensed' ]
    keys2 = ['mwarmhot', 'mhot', 'mcondensed' ]
    key_fmts2s = {}
    for k in keys2:
        key_fmts2s[k] = '%.2f,%.2f,%.2f,'

    ff =  open( 'comp.csv', 'w' )
    ff.write( '%s,%s,'%('idx_cre', 'idx') )
    for k in keys:
        ff.write( "%s,%s,%s,"%(k+'_cre', k, k+'_r') )
    for k in keys2:
        ff.write( "%s,%s,%s,"%(k+'_cre', k, k+'_r') )
    ff.write( '\n' )
    um = 1e3
    for idx_i in range(len(idxs)):
        i = idxs_cre[idx_i]
        j = idxs[idx_i]
        ff.write( '%i,%i,'%(i, j) )
        for k in keys:
            v_cre = g_cre.loc[i,k] / um
            v = g.loc[j,k] / um
            fmt = key_fmts[k]
            ff.write( fmt%(v_cre, v, (v_cre-v)/v * 100))
        for k in keys2:
            v = g.loc[j,k] / g.loc[j, 'mgas']
            v_cre = g_cre.loc[i,k] / g_cre.loc[i, 'mgas']
            fmt = key_fmts2s[k]
            ff.write( fmt%(v_cre, v, (v_cre-v)/v * 100))

        ff.write( '\n' )
    ff.close()
    
    fs = 2
    cols = 5
    rows = 6
    N = cols * rows

    dx = 1 / cols
    dy = 1 / rows

    fig = plt.figure( figsize=(fs/dx, fs/dy) )
    fig_cre = plt.figure( figsize=(fs/dx, fs/dy) )
    axs = []
    axs_cre = []
    for i in range(rows):
        for j in range(cols):
            if i * cols +  j >= len(idxs):
                break
            axs.append(fig.add_axes( [ j*dx, (rows-1-i)*dy, dx, dy ] ))
            axs_cre.append(fig_cre.add_axes( [ j*dx, (rows-1-i)*dy, dx, dy ] ))

    for i in range( N ):

        if i >= len(idxs):
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
        ax.text( 0.3*b, 0.03*a, r'$%.2f \, h^{-1}{\rm Mpc}$'%(L) )
        ax_cre.text( 0.3*b, 0.03*a, r'$%.2f \, h^{-1}{\rm Mpc}$'%(L_cre) )
    
        #mx = real2tex( m_tot, b=10, n=2 )
        #ax.text( 0.03*b, 0.9*a, r'$M: %s\, h^{-1} M_{\odot}$'%(mx) )
        #mx = real2tex( m_tot_cre, b=10, n=2 )
        #ax_cre.text( 0.03*b, 0.9*a, r'$M: %s\, h^{-1} M_{\odot}$'%(mx) )

        ax.grid()
        ax_cre.grid()

        ax.invert_yaxis()
        ax_cre.invert_yaxis()

        ax.set_xticks( [] )
        ax.set_yticks( [] )
        ax_cre.set_xticks( [] )
        ax_cre.set_yticks( [] )


    rows = len(idxs_failed) // cols
    if len(idxs_failed) % cols:
        rows += 1
    N = rows*cols
    dy = 1/rows

    fig_failed = plt.figure( figsize=(fs/dx, fs/dy) )

    axs_failed = []
    for i in range(rows):
        for j in range(cols):
            axs_failed.append(fig_failed.add_axes( [ j*dx, (rows-1-i)*dy, dx, dy ] ))

    for i in range(N):
        if i >= len(idxs_failed):
            continue
        d_failed = np.loadtxt( g_cre_dir +\
            '/Density/Density_%03i_%04i_y.dat'%(snapindex, idxs_failed[i]) )

        ax_failed = axs_failed[i]
        L_failed = np.abs(d_failed[0,1]) * 2 / 1000
        m_tot_failed= d_failed[0,8:14].sum()
        d_failed = d_failed[1:,:]
        ax_failed.imshow( d_failed, norm=mplc.LogNorm(), cmap=cm.hot )
        ax_failed.text( 0.03*b, 0.03*a, r'$L:%.2f \, h^{-1}{\rm Mpc}$'%(\
            L_failed) )
        mx = real2tex( m_tot_failed, b=10, n=2 )
        ax_failed.text( 0.03*b, 0.9*a, r'$M: %s\, h^{-1} M_{\odot}$'%(mx) )

        ax_failed.grid()
        ax_failed.invert_yaxis()
        ax_failed.set_xticks( [] )
        ax_failed.set_yticks( [] )


    fig.savefig( "g_%03i.pdf"%snapindex ) 
    fig_cre.savefig( "g_cre_%03i.pdf"%snapindex )      
    fig_failed.savefig( "g_failed_%03i.pdf"%snapindex )      


#comp()
#comp_gas()
comp_single()
