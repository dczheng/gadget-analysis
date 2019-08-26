#include "allvars.h"
#include "drfftw.h"


double *PS, binfac, Kmin, Kmax, BoxSize, *Kbin;
int bins, KN;

double powerspec_interp( double k ) {

    int ki;
    double dk;

    if ( k < Kbin[0] || k > Kbin[KN-1] )
        return 0;

    if ( k == Kbin[0] )
        //return PS[0];
        return PS[0] * CUBE( BoxSize );

    ki = 0;
    while( Kbin[ki] < k )
        ki++;
    ki--;

    dk = ( k - Kbin[ki] ) / ( Kbin[ki+1] - Kbin[ki] );

    //return ( PS[ki] * ( 1-dk ) +  PS[ki+1] * dk );
    return ( PS[ki] * ( 1-dk ) +  PS[ki+1] * dk ) * CUBE( BoxSize );

}

void test_powerspec_interp() {

    double dlogk, k;
    int N, i;
    FILE *fd;

    printf( "test powerspec_interp ...\n" );
    N = 100;
    dlogk = log( Kmax / Kmin ) / ( N-1 );
    fd = myfopen(  "w", "test_powerspec_interp.dat" );

    for ( i=0; i<N; i++ ) {

        k = i * dlogk + log( Kmin );
        k = exp( k );

        fprintf( fd, "%g %g\n", k, powerspec_interp( k ) * CUBE( BoxSize ) );

    }

    fclose( fd );
}

void compute_sigma8() {

    double R;
    sigma_struct ss;

    writelog( "compute sigma8 ...\n" );

    R = 8000;

    ss.P = &powerspec_interp;
    ss.filter = &top_hat_filter;
    ss.FilterKLimit = &top_hat_filter_k_limit;
    ss.FilterParams = &R;
    ss.norm = 1 / (sqrt(2) * M_PI);

    All.Sigma8 = sigma( ss );

    printf( "Sigma8: %g\n", All.Sigma8 );

}

void powerspec() {

    double fac, dis[6], mass, mass_tot, K[3], K2, k,
           f[3], ff, smth, po, ponorm,
           *SumPS, *SumPSUncorr, ktmp,
           *PSUncorr, *Delta, *DeltaUncorr, ComputedKmin, ComputedKmax, t, v;
    long p, index;
    int i, x, y, z, m, xyz[6], NGrid, NGrid3,
        *CountModes;;
    double *pos;
    FILE *fd;

    if ( ThisTask_Local != 0 )
        return;

    fftw_real *rhogrid;
    fftw_complex *fft_of_rhogrid;

    rfftwnd_plan fft_plan;

    writelog( "start compute power spectrum ...\n" );

    BoxSize = BoxSize;
    NGrid = All.PowSpecNGrid;
    NGrid3 = 2 * ( NGrid / 2 + 1 );
    mymalloc2( rhogrid, SQR(NGrid) * NGrid3 * sizeof(fftw_real) );
    fft_of_rhogrid = ( fftw_complex* )rhogrid;

    writelog( "particle type: ` " );
    for( i=0; i<6; i++ ) {
        if ( (1<<i) & All.PowSpecPartType )
            writelog( "%i ", i );
    }
    writelog( "` are used.\n" )

    writelog( "NGrid: %i\n", NGrid );

    mass_tot = 0;
    fac = NGrid / BoxSize;
    for( p=0; p<NumPart; p++ ) {

        if ( !( (1<<P[p].Type) & All.PowSpecPartType ) )
            continue;

        pos = P[p].Pos;

        for( i=0; i<3; i++ )
            if ( pos[i] < 0 || pos[i] > BoxSize ){
                printf( "Some Error Appearing!!!\n" );
                endrun( 20181013 );
            }

        for( i=0; i<3; i++ ) {

            t = pos[i] * fac;

            xyz[i] = t;
            xyz[i+3] = xyz[i] + 1;
            xyz[i+3] %= NGrid;

            dis[i+3] = t - xyz[i];
            dis[i] = 1 - dis[i+3];

        }

        mass_tot += P[p].Mass;
        mass = P[p].Mass;

        for( i=0; i<8; i++ ) {

            index = xyz[(i&1)*3] * NGrid * NGrid3
                  + xyz[1+((i>>1)&1)*3] * NGrid
                  + xyz[2+((i>>2)&1)*3];

            v = dis[(i&1)*3]
              * dis[1+((i>>1)&1)*3]
              * dis[2+((i>>2)&1)*3];

            rhogrid[ index ] += mass * v;
            
        }

    }

    fft_plan = rfftw3d_create_plan( NGrid, NGrid, NGrid,
            FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE );

    rfftwnd_one_real_to_complex( fft_plan, rhogrid, fft_of_rhogrid );   // in_place, so fftw_of_rhogrid is ignore.


    fac = 1.0 / mass_tot;
    Kmin = 2 * M_PI / BoxSize;
    Kmax = Kmin * BoxSize / SofteningTable[1];

    writelog( "mass tot: %g\n", mass_tot );
    writelog( "[initial] Kmin: %g, Kmax: %g [h/Mpc]\n", Kmin*1000, Kmax*1000 );

    //K1 = 2 * M_PI * BoxSize / NGrid;
    //
    bins = All.PowSpecBins;
    binfac = bins / ( log(Kmax) - log(Kmin) );
    mymalloc2( PS, sizeof( double ) * bins );
    mymalloc2( SumPS, sizeof( double ) * bins );
    mymalloc2( PSUncorr, sizeof( double ) * bins );
    mymalloc2( Delta, sizeof( double ) * bins );
    mymalloc2( SumPSUncorr, sizeof( double ) * bins );
    mymalloc2( Kbin, sizeof( double ) * bins );
    mymalloc2( Delta, sizeof( double ) * bins );
    mymalloc2( DeltaUncorr, sizeof( double ) * bins );
    mymalloc2( CountModes, sizeof( int ) * bins );

    ComputedKmin = 1e10;
    ComputedKmax = -1;
    for( x=0; x<NGrid; x++ )
        for( y=0; y<NGrid; y++ )
            for( z=0; z<NGrid; z++ ) {

                if ( x > NGrid / 2 )
                    K[0] = x - NGrid;
                else
                    K[0] = x;

                if ( y > NGrid / 2 )
                    K[1] = y - NGrid;
                else
                    K[1] = y;

                if ( z > NGrid / 2 )
                    K[2] = z - NGrid;
                else
                    K[2] = z;

                K2 = K[0] * K[0] + K[1] * K[1] + K[2] * K[2];

                if ( K2 > 0 ) {
                    if ( K2 < ( NGrid / 2.0 ) * ( NGrid / 2.0 ) ) {
                            // do deconvolution

                            f[0] = f[1] = f[2] = 1;
                            for( m=0; m<3; m++ )
                                if ( K[m] != 0 ) {
                                    f[m] = ( M_PI * K[m] ) / NGrid;
                                    f[m] = sin( f[m] ) / f[m];
                                }

                            ff = 1 / ( f[0] * f[1] * f[2] );
                            smth = ff * ff * ff * ff;
                            // end deconvolution
                        if ( z >= NGrid3/2 )
                            index = (NGrid3/2) * ( x * NGrid + y ) + ( NGrid - z );
                        else
                            index = (NGrid3/2) * ( x * NGrid + y ) + z;

                        po = ( SQR( fft_of_rhogrid[index].re ) + SQR( fft_of_rhogrid[index].im ) );

                        po *= fac * fac * smth;

                        k = sqrt( K2 ) * 2 * M_PI / BoxSize;

                        if ( k > ComputedKmax )
                            ComputedKmax = k;

                        if ( k < ComputedKmin )
                            ComputedKmin = k;

                        if ( k >= Kmin && k < Kmax ) {

                            i = log( k/Kmin ) * binfac;

                            ponorm = po / PowerSpec_Efstathiou( k );

                            if ( isnan( ponorm ) || isinf( ponorm ) ) {
                                printf( "#Error# k: %g, po: %g, Efs: %g, fft: ( %g, %g ), smth: %g, f: ( %g, %g, %g )\n",
                                        k, po, PowerSpec_Efstathiou( k ),
                                        fft_of_rhogrid[index].re,
                                        fft_of_rhogrid[index].im,
                                        smth, f[0], f[1], f[2]
                                        );
                                endrun( 20181016 );
                            }

                            //printf( "%g, %g\n", k, ponorm );
                            SumPS[i] += ponorm;
                            SumPSUncorr[i] += po;

                        CountModes[i] += 1;

                        }
                    }
                }
                /*
                printf( "( %i, %i, %i ): %g %g\n",
                            x, y, z,
                            fft_of_rhogrid[index].re,
                            fft_of_rhogrid[index].im );
                            */
            }

    writelog( "[fft] kmin: %g, kmax: %g [h/Mpc]\n", ComputedKmin*1000, ComputedKmax*1000 );

    index = 0;

    ktmp = Kmin;

    Kmin = 1e10;
    Kmax= -1;

    for( i=0; i<bins; i++ ) {

        k = exp( (i+0.5) / binfac + log( ktmp ) );

        if ( CountModes[i] > 0 ) {

            PS[index] = PowerSpec_Efstathiou( k ) * SumPS[i] / CountModes[i];
            PSUncorr[index] = SumPSUncorr[i] / CountModes[i];

            Kbin[index] = k;

            if ( k < Kmin )
                Kmin = k;

            if ( k > Kmax )
                Kmax = k;

            CountModes[index] = CountModes[i];
            Delta[index] = 4 * M_PI * pow( k, 3 ) / pow( 2 * M_PI / BoxSize, 3 ) * PS[index];
            DeltaUncorr[index] = 4 * M_PI * pow( k, 3 ) / pow( 2 * M_PI / BoxSize, 3 ) * PSUncorr[index];
            index ++;
        }
    }

    KN = index;

    writelog( "[final] Kmin: %g, Kmax: %g [h/Mpc]\n", Kmin*1000, Kmax*1000 );

    create_dir( "./PowSpec" );
    fd = myfopen( "w", "./PowSpec/PowSpec_%03i.dat", SnapIndex );

    for( i=0; i<KN; i++ ) {
        fprintf( fd, "%g %g %g %g %g %g\n",
            Kbin[i],
            Delta[i], DeltaUncorr[i],
            (double)CountModes[i],
            PS[i], PSUncorr[i] );
    }

    fclose( fd );

    rfftwnd_destroy_plan( fft_plan );

    //test_powerspec_interp();
    //
    if ( Kmin*1000 > 1e-2 )
        writelog( "Kmin is too small, sigma8 can't be computed!\n" )
    else
        compute_sigma8();

    myfree( PS );
    myfree( SumPS );
    myfree( PSUncorr );
    myfree( SumPSUncorr );
    myfree( Delta );
    myfree( CountModes );
    myfree( Kbin );
    myfree( Delta );
    myfree( DeltaUncorr );
    myfree( rhogrid );

    writelog( "compute power spectrum ... done.\n" );
    put_sep;

}
