#include "allvars.h"

double kernel_func( x ) {
    return ( 1 / sqrt( 2 * PI ) / All.ConvSigma *
        exp( -0.5 * SQR(  x / All.ConvSigma ) ) );
}

void init_conv_kernel( double ds ) {
    int i, j, c, N;
    double r;
    writelog( "init conv kernel ...\n" );
    if ( All.ConvN % 2 == 0 )
        endruns( "Even ConvN is not supported!" );

    mymalloc1( ConvKernel, SQR( All.ConvN ) * sizeof( double ) );
    N = All.ConvN;
    c = N / 2 + 1;
    for ( i=0; i<N; i++ )
        for ( j=0; j<N; j++ ) {
            r = sqrt( SQR( i-c ) + SQR( j-c ) ) * ds;
            ConvKernel[ i*N + j ] = kernel_func( r );
        }
    writelog( "init conv kernel ... done.\n" );
}

void free_conv_kernel() {
    writelog( "free conv kernel ....\n" );
    myfree( ConvKernel );
}

void conv( double ds ) {
    int i, j, ii, jj, N,  k, l;
    double *img;
    writelog( "conv ...\n" );
    N = All.ConvN;
    mymalloc2( img, PicSize2 * sizeof( double ) );
    for ( i=0; i<PicSize; i++ )
        for ( j=0; j<PicSize; j++ )
            for ( ii=0; ii<N; ii++ )
                for ( jj=0; jj<N; jj++ ) {
                    k = i+ii-N/2+1;
                    l = j+jj-N/2+1;
                    if ( k < 0 ||
                         k > PicSize ||
                         l< 0 ||
                         l> PicSize )
                        continue;
                    img[ i*PicSize + j ] += ConvKernel[ ii*N + jj ] *
                         image.img[ k * PicSize + l ];
                }
    memcpy( image.img, img, PicSize2 * sizeof( double ) );
    myfree( img );
    writelog( "conv ... done.\n" );
}

