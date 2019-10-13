#include "allvars.h"
#include "drfftw.h"

void pm_potential( const double *data, long N, double L, double NGrid, int mode, double *pot ) {

    // mode 0: non-periodic  mode 1: periodic

    int i, j, k, m, NGrid3, ks[3], k2;
    long ip, ip1;
    double fac, fac_pot, f[3], smth, ff, x[3], r, u, re, im;

    fftw_real *rhogrid, *kernel;
    fftw_complex *fft_of_rhogrid, *fft_of_kernel;
    rfftwnd_plan fft_plan, fft_plan_inv;
    NGrid3 = 2 * ( NGrid / 2 + 1  );

    mymalloc2( rhogrid, SQR(NGrid) * NGrid3 * sizeof(fftw_real)  );
    fft_of_rhogrid = ( fftw_complex*  )rhogrid;

    fft_plan = rfftw3d_create_plan( NGrid, NGrid, NGrid,
            FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE ); 
    fft_plan_inv = rfftw3d_create_plan( NGrid, NGrid, NGrid,
            FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE ); 


    if ( !mode  ) {
        mymalloc2( kernel, SQR(NGrid) * NGrid3 * sizeof(fftw_real)  );
        fft_of_kernel = ( fftw_complex*  )kernel;
        fac_pot = G / pow( L, 4 ) * pow( L/NGrid, 3 );
    }
    else {

        fac_pot = G / ( PI * L );
        /*
        to get potential
        [4 * PI * G / ( k_i * 2*Pi / L )^2]  * [1 / ( L/N )^3] * N^3,
        1 / ( L/N )^3: mass to density;
        N^3: implementation of fftw, that is FFTW computes an unnormalized transform,
        */
        /*
        asmth = ASMTH * ( L/NGrid );
        asmth2 = ( 2*PI ) / L * asmth;
        asmth2 *= asmth2;
        */

    }

    if ( !mode ) {
        for( i=0; i<NGrid; i++ )
            for ( j=0; j<NGrid; j++ )
                for( k=0; k<NGrid; k++ ) {
                    x[0] = ((double)i) / NGrid;
                    x[1] = ((double)j) / NGrid;
                    x[2] = ((double)k) / NGrid;
    
                    for( m=0, r=0; m<3; m++ ) {
                        x[m] = ( x[m] >= 0.5 ? x[m]-1.0 : x[m] );
                        r += x[m] * x[m];
                    }
                    r = sqrt(r);
                    u = 0.5 * r / ( ((double)ASMTH) / NGrid );
                    fac = 1 - erfc(u);
                    ip = i * NGrid * NGrid3 + j * NGrid3 + k;
                    if ( r > 0 )
                        kernel[ip] = -fac/r;
                    else
                        kernel[ip] = -1 / ( sqrt(PI) * (((double)ASMTH)/NGrid));
    
                }
    
        rfftwnd_one_real_to_complex( fft_plan, kernel, NULL  );

        for( i=0; i<NGrid; i++ )
            for( j=0; j<NGrid; j++ )
                for( k=0; k<NGrid/2+1; k++ ) {
    
                    if ( i==0 && j==0 && k==0 ) {
                        continue;
                    }
    
                    ks[0] = ( i>NGrid/2 ) ? i-NGrid : i;
                    ks[1] = ( j>NGrid/2 ) ? j-NGrid : j;
                    ks[2] = ( k>NGrid/2 ) ? k-NGrid : k;
    
                    for( m=0, k2=0; m<3; m++ ) {
                        k2 += ks[m]*ks[m];
                        f[m] = 1;
                        if ( ks[m] != 0 ) {
                            f[m] = ( PI * ks[m] ) / NGrid;
                            f[m] = sin(f[m]) / f[m];
                        }
                    }
    
                    ff = 1 / ( f[0] * f[1] * f[2] );
                    ff = ff * ff * ff * ff;
    
                    ip = i * NGrid * ( NGrid/2+1 ) + j * ( NGrid/2+1 ) + k;
                    fft_of_kernel[ ip ].re *= ff; 
                    fft_of_kernel[ ip ].im *= ff;
                }
    }

    
    field_to_grid_fft( data, rhogrid, L, NGrid, N, 1, mode );

    rfftwnd_one_real_to_complex( fft_plan, rhogrid, NULL  );

    for( i=0; i<NGrid; i++ )
        for( j=0; j<NGrid; j++ )
            for( k=0; k<NGrid/2+1; k++ ) {

                if ( !mode ) {
                    ip = i * NGrid * ( NGrid/2+1 ) + j * ( NGrid/2+1 ) + k;
                    re = fft_of_rhogrid[ip].re * fft_of_kernel[ip].re -
                        fft_of_rhogrid[ip].im * fft_of_kernel[ip].im;

                    im = fft_of_rhogrid[ip].re * fft_of_kernel[ip].im +
                        fft_of_rhogrid[ip].im * fft_of_kernel[ip].re;

                    fft_of_rhogrid[ip].re = fac_pot * re;
                    fft_of_rhogrid[ip].im = fac_pot * im;
                    continue;
                }

                /* periodic */
                if ( i==0 && j==0 && k==0 ) {
                    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0;
                    continue;
                }

                ks[0] = ( i>NGrid/2 ) ? i-NGrid : i;
                ks[1] = ( j>NGrid/2 ) ? j-NGrid : j;
                ks[2] = ( k>NGrid/2 ) ? k-NGrid : k;

                for( m=0, k2=0; m<3; m++ ) {
                    k2 += ks[m]*ks[m];
                    f[m] = 1;
                    if ( ks[m] != 0 ) {
                        f[m] = ( PI * ks[m] ) / NGrid;
                        f[m] = sin(f[m]) / f[m];
                    }
                }

                ff = 1 / ( f[0] * f[1] * f[2] );
                //smth = - exp( -k2 * asmth2 ) / k2 * fac_pot * ff * ff * ff * ff;
                smth = - 1 / k2 * fac_pot * ff * ff * ff * ff;

                ip = i * NGrid * ( NGrid/2+1 ) + j * ( NGrid/2+1 ) + k;
                fft_of_rhogrid[ ip ].re *= smth; 
                fft_of_rhogrid[ ip ].im *= smth;

            }

    rfftwnd_one_complex_to_real( fft_plan_inv, fft_of_rhogrid, NULL  );

    for( i=0; i<NGrid; i++ )
        for ( j=0; j<NGrid; j++ )
            for( k=0; k<NGrid; k++ ) {
                ip  = i * NGrid * NGrid3 + j * NGrid3 + k;
                ip1 = i * NGrid * NGrid + j * NGrid + k;
                pot[ip1] = rhogrid[ip];
            }

    rfftwnd_destroy_plan( fft_plan  );
    rfftwnd_destroy_plan( fft_plan_inv );
    myfree( rhogrid );

    if ( !mode )
        myfree( kernel );

}
