#!/usr/bin/env python3

from my_work_env import *
import matplotlib.gridspec as gridspec
import matplotlib.ticker as tik

def plot_pdf1():
    fn = "cre_256_no_turb_CrePressurePdf_0.00.dat"

    ds = np.loadtxt( fn )
    data_info = ds[ 0, : ]
    ds = ds[1:,:]

    XMin = data_info[1]
    XMax = data_info[2]
    YMin = data_info[3]
    YMax = data_info[4]
    z = data_info[10]


    m,n = ds.shape
    xmin = int( XMin ) + 1
    xmax = int( XMax )
    xl = XMax - XMin
    ymin = int( YMin ) + 1
    ymax = int( YMax )
    yl = YMax - YMin

    print( "xmin: %i, xmax: %i, ymin: %i, ymax: %i"%( xmin, xmax, ymin ,ymax ) )

    #xLocList = np.linspace( xmin, xmax, xmax - xmin + 1 )
    xLocList = np.arange( xmin, xmax+1, 2 )
    xFmtList = [ r"$10^{%.0f}$"%i for i in xLocList ]
    xLocList = ( xLocList - XMin ) / xl * (n-1)
    xLocList = np.hstack( [xLocList, n+xLocList, 2*n+xLocList] )
    print( "xticks: ", xFmtList )

    #yLocList = np.linspace( ymin, ymax, ymax - ymin + 1 )
    yLocList = np.arange( ymin, ymax+1, 2 )
    yFmtList = [ r"$10^{%.0f}$"%i for i in yLocList ]
    yLocList = ( yLocList - YMin ) / yl * (m-1)
    print( "yticks: ", yFmtList )

    xloc = tik.FixedLocator(xLocList)
    xfmt = tik.FixedFormatter(xFmtList)
    yloc = tik.FixedLocator(yLocList)
    yfmt = tik.FixedFormatter(yFmtList)

    my_norm = mplc.LogNorm()
    my_cmap = cm.jet

    img = plt.imshow( ds, norm=my_norm, cmap=my_cmap )
    #cbar = plt.colorbar( img )
    cbar = plt.colorbar( img, pad=0)
    cbar.set_label( r'$f(R_{\rm CR}, R_{\rm CRE})$', fontsize=20 )
    #cbar.ax.tick_params( direction='in', pad=5, labelsize=20 )
    cbar.ax.tick_params( direction='in', width=0.3, length=1.5, labelsize=20 )

    #t = cbar.ax.get_yticklabels()
    #cbar.ax.set_yticklabels( t, fontsize=7 )
        #exit()

    ax = plt.gca()
    fig = plt.gcf()

    ax.invert_yaxis()

    ax.xaxis.set_major_locator( xloc )
    ax.xaxis.set_major_formatter( xfmt )
    ax.xaxis.set_tick_params(labelsize=5)

    ax.yaxis.set_major_locator( yloc )
    ax.yaxis.set_major_formatter( yfmt )
    ax.yaxis.set_tick_params(labelsize=5)

    ax.tick_params( axis='both', direction='in', pad=5, labelsize=20 )

    ax.set_ylabel( r"$R_{\rm CRE}$", fontsize=20 )
    ax.set_xlabel( r"$R_{\rm CR}$", fontsize=20 )
    #ax.set_title( "CRE Pressure PDF" )

    plt.tight_layout()

    plt.savefig( output_dir + 'PPdf.pdf', figsize=(5,5) )

def gen_bin( x, y ):
    xx = []
    yy = []
    y = y / y.sum() * 100
    for i in range( len(x)-1 ):
        xx.append( x[i] )
        xx.append( x[i+1] )
        yy.append( y[i] )
        yy.append( y[i] )

    return ( xx, yy )

def plot_pdf2():
    fn = sys.argv[1]
    fn_out = sys.argv[2]

    d = np.loadtxt( fn );

    fig, ax = plt.subplots( 1, 1 )

    x, y = gen_bin( d[:,0], d[:,1] )
    ax.plot( x, y, label=r'$\frac{P_{\rm CRE}}{P_{\rm tot}}$' )

    x, y = gen_bin( d[:,2], d[:,3] )
    ax.plot( x, y, label=r'$\frac{P_{\rm CRP}}{P_{\rm tot}}$' )

    #ax.set_yscale( 'log' )
    ax.set_xscale( 'log' )
    ax.set_ylabel( r'$\rm Percent$', fontsize=20 )
    ax.tick_params( axis='both', direction='in', pad=5, labelsize=10 )
    ax.legend( prop={'size': 20}, framealpha=0.1 )
    ax.minorticks_off()
    ax.locator_params(numticks=12)

    ylim = list( ax.get_ylim() )
    ylim[0] = 1e-1
    ax.set_ylim( ylim )

    xlim = list( ax.get_xlim() )
    xlim[0] = 1e-12
    ax.set_xlim( xlim )

    plt.savefig( fn_out )

plot_pdf2()
