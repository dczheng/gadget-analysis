#include "allvars.h"

double kernel( double q ) {
    if ( 0 <= q && q <= 0.5 )
        return 8.0 / PI * ( 1 - 6*q*q + 6*q*q*q );
    if ( q <= 1 )
        return 8.0 / PI * ( 2.0 * pow( 1-q, 3.0 ) );
    return 0;
}

void init_kernel_matrix() {
    int pt, i, j, k, N, Nhalf;
    long num;
    double dh, q, h;
    writelog( "initialize kernel matrix ...\n" );
    N = All.KernelN;
    Nhalf = N / 2;
    for( pt=0; pt<6; pt++ ) {
        num = header.npartTotal[pt] + ( ( ( long )header.npartTotalHighWord[pt] ) << 32 );
        if ( num != 0 ) {
            KernelMat2D[pt] = malloc( sizeof( double ) * N * N );
            KernelMat3D[pt] = malloc( sizeof( double ) * N * N * N );
            memset( KernelMat2D[pt], 0, sizeof( double ) * N * N );
            memset( KernelMat3D[pt], 0, sizeof( double ) * N * N * N );

            if ( SofteningTable[pt] == 0 )
                endrun0( "SofteningTable[%d] is zeros !!!\n", pt );

            h = SofteningTable[pt];
            dh = h / Nhalf;
            for ( i=0; i<N; i++ )
                for ( j=0; j<N; j++ )
                    for ( k=0; k<N; k++ ) {
                        //writelog( "%i\n", Nhalf );
                        q = sqrt( SQR( i-Nhalf ) +
                                  SQR( j-Nhalf ) +
                                  SQR( k-Nhalf ) );
                        q = q / Nhalf;
                        KernelMat3D[pt][ i*N*N + j*N + k ] = kernel( q ) / pow( h, 3 )
                            * pow( dh, 3 );
                        KernelMat2D[pt][ i*N + j ] += KernelMat3D[pt][ i*N*N + j*N + k ];
                        //writelog( "%f, %f\n", q, KernelMat2D[pt][i*N+j] );
                    }
            /*
            double t;
            t = 0;
            for ( i=0; i<N; i++ )
                for ( j=0; j<N; j++ ) {
                    t += KernelMat2D[pt][ i*N + j ];
                }
            writelog( "%f\n", t );
            endrun( 0 );
            */

            for ( i=0, q=0; i<N*N*N; i++ ) {
                q += KernelMat3D[pt][i];
            }
            writelog( "sum of 3d kernel[%d] = %g\n", pt, q );
            for ( i=0, q=0; i<N*N; i++ )
                q += KernelMat2D[pt][i];
            writelog( "sum of 2d kernel[%d] = %g\n", pt, q );
        }
    }
    writelog( "initialize kernel matrix ... done.\n" );
    put_sep0;
}

void free_kernel_matrix() {
    long num;
    int pt;
    writelog( "free kernel matrix ...\n" );
    for ( pt=0; pt<6; pt++ ) {
        num = header.npartTotal[pt] + ( ( ( long )header.npartTotalHighWord[pt] ) << 32 );
        if ( num != 0 ) {
            free( KernelMat2D[pt] );
            free( KernelMat3D[pt] );
        }
    }
}

