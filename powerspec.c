#include "allvars.h"
#include "srfftw.h"

void powerspec() {

    double K0, K1, fac, dis[3], mass, mass_tot, K[3], K2, k, BoxSize,
           f[3], ff, smth, po, ponorm, binfac, *SumPS, *PS, *SumPSUncorr,
           *PSUncorr, *Kbin, *Delta, *DeltaUncorr, kmin, kmax;
    long p;
    int i, x, y, z, m, xyz[3], xyz1[3], NGrid, NGrid3, index,
        bins, *CountModes;;
    MyFloat *pos;
    FILE *fd;
    char buf[100];

    fftw_real *rhogrid;
    fftw_complex *fft_of_rhogrid;


    rfftwnd_plan fft_plan;

    writelog( "start compute power spectrum ...\n" );

    BoxSize = All.BoxSize;
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

            /*
            while( pos[i] < 0 )
                pos[i] += BoxSize;

            while( pos[i] > BoxSize )
                pos[i] -= BoxSize;
                */

            dis[i] = pos[i] * fac;
            xyz[i] = dis[i];
            dis[i] -= xyz[i];

            //printf( "%i, dis: %g, xyz: %i\n", i, dis[i], xyz[i] );
            xyz1[i] = xyz[i] + 1;

            xyz1[i] %= NGrid;

        }

        //printf( "( %i, %i, %i )\n", xyz[0], xyz[1], xyz[2] );

        mass_tot += P[p].Mass;
        mass = P[p].Mass;

        if ( P[p].Type == 0 ) {

            if ( 1<<4 & All.PowSpecPartType ) {
                mass_tot += SphP[p].Star_Mass;
                mass += SphP[p].Star_Mass;
            }

            if ( 1<<5 & All.PowSpecPartType ) {
                mass += SphP[p].BH_Mass;
                mass_tot += SphP[p].BH_Mass;
            }

        }

//        rhogrid[ NGrid3 * ( NGrid * xyz[0] + xyz[1] ) + xyz[2] ] += mass;
        rhogrid[ NGrid3 * ( NGrid * xyz[0]  + xyz[1]  ) + xyz[2]  ]  += mass * ( 1-dis[0] ) * ( 1-dis[1] ) * ( 1-dis[2] );
        rhogrid[ NGrid3 * ( NGrid * xyz[0]  + xyz[1]  ) + xyz1[2] ]  += mass * ( 1-dis[0] ) * ( 1-dis[1] ) * dis[2];
        rhogrid[ NGrid3 * ( NGrid * xyz[0]  + xyz1[1] ) + xyz[2]  ]  += mass * ( 1-dis[0] ) * dis[1]       * ( 1-dis[2] );
        rhogrid[ NGrid3 * ( NGrid * xyz[0]  + xyz1[1] ) + xyz1[2] ]  += mass * ( 1-dis[0] ) * dis[1]       * dis[2];
        rhogrid[ NGrid3 * ( NGrid * xyz1[0] + xyz[1]  ) + xyz[2]  ]  += mass * dis[0]       * ( 1-dis[1] ) * ( 1-dis[2] );
        rhogrid[ NGrid3 * ( NGrid * xyz1[0] + xyz[1]  ) + xyz1[2] ]  += mass * dis[0]       * ( 1-dis[1] ) * dis[2];
        rhogrid[ NGrid3 * ( NGrid * xyz1[0] + xyz1[1] ) + xyz[2]  ]  += mass * dis[0]       * dis[1]       * ( 1-dis[2] );
        rhogrid[ NGrid3 * ( NGrid * xyz1[0] + xyz1[1] ) + xyz1[2] ]  += mass * dis[0]       * dis[1]       * dis[2];

    }

    /*
    for( i=0; i<dim[0]; i++ )
        for( j=0; j<NGrid; j++ )
            for( k=0; k<NGrid3; k++ )
                rhogrid[ NGrid3 * ( NGrid * i + j ) + 0 ] += rhogrid[ NGrid3 * ( NGrid * i + j ) + k ];

    fd = fopen( "rhogrid.dat", "w" );
    for( i=0; i<dim[0]; i++ ) {
        for( j=0; j<NGrid; j++ )
            fprintf( fd, "%g ", rhogrid[ NGrid3 * ( NGrid * i + j ) ] );
        fprintf( fd, "\n" );
    }
    fclose( fd );
    return;
    */

    fft_plan = rfftw3d_create_plan( NGrid, NGrid, NGrid,
            FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE );

    rfftwnd_one_real_to_complex( fft_plan, rhogrid, fft_of_rhogrid );   // in_place, so fftw_of_rhogrid is ignore.


    fac = 1.0 / mass_tot;
    K0 = 2 * M_PI / BoxSize;
    K1 = K0 * BoxSize / All.SofteningTable[1];

    writelog( "mass tot: %g\n", mass_tot );
    writelog( "K0: %g, K1: %g\n", K0, K1 );
    //K1 = 2 * M_PI * BoxSize / NGrid;
    //
    bins = All.PowSpecBins;
    binfac = bins / ( log(K1) - log(K0) );
    mymalloc2( PS, sizeof( double ) * bins );
    mymalloc2( SumPS, sizeof( double ) * bins );
    mymalloc2( PSUncorr, sizeof( double ) * bins );
    mymalloc2( Delta, sizeof( double ) * bins );
    mymalloc2( SumPSUncorr, sizeof( double ) * bins );
    mymalloc2( Kbin, sizeof( double ) * bins );
    mymalloc2( Delta, sizeof( double ) * bins );
    mymalloc2( DeltaUncorr, sizeof( double ) * bins );
    mymalloc2( CountModes, sizeof( int ) * bins );

    kmin = 1e10;
    kmax = -1;
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

                        if ( k > kmax )
                            kmax = k;

                        if ( k < kmin )
                            kmin = k;

                        if ( k >= K0 && k < K1 ) {

                            i = log( k/K0 ) * binfac;

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

    writelog( "[fft] kmin: %g, kmax: %g\n", kmin, kmax );
    for( i=0; i<bins; i++ ) {

        Kbin[i] = exp( (i+0.5) / binfac + log( K0 ) );

        if ( CountModes[i] > 0 ) {
            PS[i] = PowerSpec_Efstathiou( Kbin[i] ) * SumPS[i] / CountModes[i];
            PSUncorr[i] = SumPSUncorr[i] / CountModes[i];
        }
        else {
            PS[i] = 0;
            PSUncorr[i] = 0;
        }

        Delta[i] = 4 * M_PI * pow( Kbin[i], 3 ) / pow( 2 * M_PI / BoxSize, 3 ) * PS[i];
        DeltaUncorr[i] = 4 * M_PI * pow( Kbin[i], 3 ) / pow( 2 * M_PI / BoxSize, 3 ) * PSUncorr[i];

        /*
        printf( "[%i] %g %g %g %g %g %g\n", i, (double)CountModes[i], SumPS[i], PS[i], SumPSUncorr[i],
                Delta[i], DeltaUncorr[i] );
                */

    }

    create_dir( "./PowSpec" );
    sprintf( buf, "./PowSpec/PowSpec_%.2f.dat", All.RedShift );

    fd = fopen( buf, "w" );

    for( i=0; i<bins; i++ ) {
        fprintf( fd, "%g %g %g %g %g %g\n",
                Kbin[i],
                Delta[i], DeltaUncorr[i],
                (double)CountModes[i],
                PS[i], PSUncorr[i] );
    }

    fclose( fd );

    rfftwnd_destroy_plan( fft_plan );
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
    put_block_line;

}
