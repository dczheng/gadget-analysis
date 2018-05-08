#include "allvars.h"

double kernel_func( x ) {
    return ( 1 / sqrt( 2 * PI ) / All.ConvSigma *
        exp( -0.5 * SQR(  x / All.ConvSigma ) ) );
}

void init_conv_kernel( double ds ) {
    int i, j, c, N;
    double r;
    writelog( "init conv kernel ...\n" );
    if ( All.ConvN % 2 == 0 ){
        printf( "Even ConvN is not supported!\n" );
        endrun( 20180508 );
    }
    mymalloc( All.ConvKernel, SQR( All.ConvN ) * sizeof( double ) );
    N = All.ConvN;
    c = N / 2 + 1;
    for ( i=0; i<N; i++ )
        for ( j=0; j<N; j++ ) {
            r = sqrt( SQR( i-c ) + SQR( j-c ) ) * ds;
            All.ConvKernel[ i*N + j ] = kernel_func( r );
        }
    writelog( "init conv kernel ... done.\n" );
}

void free_conv_kernel() {
    writelog( "free conv kernel ....\n" );
    myfree( All.ConvKernel );
}
void conv( double ds ) {
    int PicSize, i, j, ii, jj, N, PicSize2, k, l;
    double *img;
    PicSize = All.PicSize;
    PicSize2 = All.PicSize2;
    writelog( "conv ...\n" );
    N = All.ConvN;
    mymalloc( img, PicSize2 * sizeof( double ) );
    memset( img, 0, PicSize2 * sizeof( double ) );
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
                    img[ i*PicSize + j ] += All.ConvKernel[ ii*N + jj ] *
                         image.img[ k * PicSize + l ];
                }
    memcpy( image.img, img, PicSize2 * sizeof( double ) );
    myfree( img );
    writelog( "conv ... done.\n" );
}

