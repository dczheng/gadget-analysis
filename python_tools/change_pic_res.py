#!/usr/bin/env python3
import numpy as np

def get_index( i, t ):
    ii = i * t
    t = int( ii )
    di = ii - t
    ii = t
    return ( ii, di )

def get_index2( i, j, t ):
    ii, di = get_index( i, t )
    jj, dj = get_index( j, t )
    ii1, di1 = get_index( i+1, t )
    jj1, dj1 = get_index( j+1, t )
    return ( [ii, jj, ii1, jj1], [di, dj, di1, dj1] )

def check_index( i, M, sig ):
    if i >=M:
        print( "[%i] bug???"%sig )
        i = M-1
    return i

def check_index2( ii, M, N, sig ):
    t = [ M, N, M, N ]
    for i in range(4):
        ii[0] = check_index( ii[0], t[i], sig*10+i )
    return ii

def get_r( d, r, flag ):
    if flag == 1:
        t1 = ( 1-d[1] ) * r
        t2 = d[3] * r
        return ( t1, t2 )

    if flag == 2:
        t1 = ( 1-d[0] ) * r
        t2 = d[2] * r
        return ( t1, t2 )

    if flag == 3:
        t1 = (1-d[0]) * (1-d[1])
        t2 =  d[2] * d[3]
        t3 = (1-d[0]) * d[3]
        t4 = d[2] * (1-d[1])
        return ( t1, t2, t3, t4 )

def gen1( d1, res1, res2 ):

    M1, N1 = d1.shape
    M2 = int(M1 * res1 / res2) - 2
    N2 = int(N1 * res1 / res2) - 2
    d2 = np.zeros( ( M2, N2 ) )

    print( "res: %.2f -> %.2f"%(res1, res2) )
    print( "(%i, %i) -> (%i,%i)"%(M1, N1, M2, N2) )

    t21 = res2 / res1
    r21 = res2 * res2 / ( res1 * res1 )

    sig = M2 * N2 // 10
    for i in range( M2 ):
        for j in range( N2 ):

            if (i*N2+j) % sig == 0:
                #print( i*M2+j )
                print( "%3.0f%%"%( (i*N2+j) / (M2*N2) * 100 ) )

            iii, dii = get_index2( i+1, j+1, t21 )
            iii = check_index2( iii, M1, N1, 1 )

            ii = iii[0]
            jj = iii[1]
            ii1 = iii[2]
            jj1 = iii[3]

            if (ii == ii1 and jj == jj1):
                d2[i,j] +=  d1[ii,jj] * r21
                continue


            if ii == ii1:
                t1, t2 = get_r( dii, t21, 1 )
                d2[i,j] += d1[ii,jj] * t1 + d1[ii,jj1] * t2
                #print( t1+t2 )
                continue

            if jj == jj1:
                t1, t2 = get_r( dii, t21, 2 )
                d2[i,j] += d1[ii,jj] * t1 + d1[ii1,jj] * t2
                #print( t1+t2 )
                continue


            t1, t2, t3, t4 = get_r( dii, t21, 3 )
            #print( t1+t2+t3+t4 )
            d2[i,j] += d1[ii,jj]   * t1
            d2[i,j] += d1[ii1,jj1] * t2
            d2[i,j] += d1[ii,jj1]  * t3
            d2[i,j] += d1[ii1,jj]  * t4

    return d2

def gen2( d1, res1, res2 ):

    M1, N1 = d1.shape
    M2 = int(M1 * res1 / res2) + 2
    N2 = int(N1 * res1 / res2) + 2
    d2 = np.zeros( ( M2, N2 ) )

    print( "res: %.2f -> %.2f"%(res1, res2) )
    print( "(%i, %i) -> (%i,%i)"%(M1, N1, M2, N2) )

    t12 = res1 / res2
    r21 = res2 * res2 / ( res1 * res1 )

    sig = M1 * N1 // 10
    for i in range( M1 ):
        for j in range( N1 ):

            if (i*N1+j) % sig == 0:
                #print( i*M2+j )
                print( "%3.0f%%"%( (i*N1+j) / (M1*N1) * 100 ) )

            iii, dii = get_index2( i, j, t12 )
            iii = check_index2( iii, M2, N2, 2 )

            ii  = iii[0] + 1
            jj  = iii[1] + 1
            ii1 = iii[2] + 1
            jj1 = iii[3] + 1

            if ii == ii1 and jj == jj1:
                d2[ii,jj] +=  d1[i,j]
                continue

            if ii == ii1:
                t1, t2 = get_r( dii, t12, 1 )
                d2[ii,jj]  += d1[i,j] * t1 * r21
                d2[ii,jj1] += d1[i,j] * t2 * r21
                continue


            if jj == jj1:
                t1, t2 = get_r( dii, t12, 2 )
                d2[ii,jj] += d1[i,j]  * t1 * r21
                d2[ii1,jj] += d1[i,j] * t2 * r21
                continue

            t1, t2, t3, t4 = get_r( dii, t12*r21, 3 )
            #print( t1+t2+t3+t4 )
            d2[ii,jj]   += d1[i,j] * t1 * r21
            d2[ii1,jj1] += d1[i,j] * t2 * r21
            d2[ii,jj1]  += d1[i,j] * t3 * r21
            d2[ii1,jj]  += d1[i,j] * t4 * r21
        continue
        a = 4
        if i >= a:
            a, dii = get_index2( a+10, j, t12 )
            for j in range( M2 ):
                if j % 2 == 0:
                    d2[a[0]:,j] = 2
                else:
                    d2[a[0]:,j] = -2
            break

    return d2

def gen( d, res1, res2 ):
    if res1 > res2:
        return gen1( d, res1, res2 )
    if res1 < res2:
        return gen2( d, res1, res2 )

    print( "res1 == res2 !!!!" )
    exit()


def test():

    import matplotlib
    matplotlib.use( 'agg' )
    import matplotlib.pyplot as plt
    #import astropy.io.fits as fits
    import matplotlib.colors as mplc
    from matplotlib import cm
    f = fits.getdata('/home/xkchen/mywork/ICL/data/test_data/frame-r-ra36.455-dec-5.896-redshift0.233.fits',header=True)
    d1 = f[0]
    #M = N = 100
    #a = 10
    #M = N = 100
    #a = 10
    #d1 = np.zeros( [M, N] )
    #d1[M//a:M-M//a,N//a:N-N//a] = np.random.rand( M-M*2//a,N-N*2//a )
    #d1[M//a:M-M//a,N//a:N-N//a] = 1

    b = 2.31
    d2 = gen( d1, b, 1 )
    d3 = gen( d1, 1, b )
    #d3 =  d2

    fig = plt.figure()
    ax1 = fig.add_subplot( 221 )
    ax2 = fig.add_subplot( 222 )
    ax3 = fig.add_subplot( 223 )
    #ax1.imshow( d1, norm=mplc.LogNorm(), cmap=cm.jet )
    #ax2.imshow( d2, norm=mplc.LogNorm(), cmap=cm.jet )
    #ax3.imshow( d3, norm=mplc.LogNorm(), cmap=cm.jet )
    cmap = cm.jet
    img = ax1.imshow( d1, cmap=cmap )
    plt.colorbar( img, ax=ax1 )

    img = ax2.imshow( d2, cmap=cmap )
    plt.colorbar( img, ax=ax2 )

    img = ax3.imshow( d3, cmap=cmap )
    plt.colorbar( img, ax=ax3 )

    #plt.show()
    fig.tight_layout()
    fig.savefig( 't.png' )

def main():
    test()

if __name__ == '__main__':
    main()
