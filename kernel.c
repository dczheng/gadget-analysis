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
            All.KernelMat2D[pt] = malloc( sizeof( double ) * N * N );
            All.KernelMat3D[pt] = malloc( sizeof( double ) * N * N * N );
            memset( All.KernelMat2D[pt], 0, sizeof( double ) * N * N );
            memset( All.KernelMat3D[pt], 0, sizeof( double ) * N * N * N );
            if ( All.SofteningTable[pt] == 0 ) {
                printf( "SofteningTable[%d] is zeros !!!\n", pt );
                endrun( 20171207 );
            }
            h = All.SofteningTable[pt];
            dh = h / Nhalf;
            for ( i=0; i<N; i++ )
                for ( j=0; j<N; j++ )
                    for ( k=0; k<N; k++ ) {
                        q = sqrt( pow( i-Nhalf, 2 ) +
                                  pow( j-Nhalf, 2 ) +
                                  pow( k-Nhalf, 2 ) );
                        q = q / Nhalf;
                        All.KernelMat3D[pt][ i*N*N + j*N + k ] = kernel( q ) / pow( h, 3 )
                            * pow( dh, 3 );
                        All.KernelMat2D[pt][ i*N + j ] += All.KernelMat3D[pt][ i*N*N + j*N + k ];
                    }
            for ( i=0, q=0; i<N*N*N; i++ ) {
                q += All.KernelMat3D[pt][i];
            }
            writelog( "sum of 3d kernel[%d] = %g\n", pt, q );
            for ( i=0, q=0; i<N*N; i++ )
                q += All.KernelMat2D[pt][i];
            writelog( "sum of 2d kernel[%d] = %g\n", pt, q );
        }
    }
    writelog( "initialize kernel matrix ... done.\n" );
}

void free_kernel_matrix() {
    long num;
    int pt;
    writelog( "free kernel matrix ...\n" );
    for ( pt=0; pt<6; pt++ ) {
        num = header.npartTotal[pt] + ( ( ( long )header.npartTotalHighWord[pt] ) << 32 );
        if ( num != 0 ) {
            free( All.KernelMat2D[pt] );
            free( All.KernelMat3D[pt] );
        }
    }
}

