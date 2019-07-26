#include "allvars.h"
#include "srfftw.h"

void corr( double *x, double *y, double *c, int N ) {

    int N3;
    long i, j, k, index, index0;
    double a, b;

    fftw_real *xx, *yy;
    fftw_complex *kxx, *kyy;
    rfftwnd_plan fft_plan, fft_plan_inv;
    N3 = 2 * ( N/2+1 );

    mymalloc2( xx, SQR(N)*N3*sizeof( fftw_real ) );
    mymalloc2( yy, SQR(N)*N3*sizeof( fftw_real ) );
    kxx = ( fftw_complex* ) xx;
    kyy = ( fftw_complex* ) yy;

    for( i=0; i<N; i++ )
        for( j=0; j<N; j++ )
                for( k=0; k<N; k++ ) {
                    index = i*N*N3 + j*N3 + k;
                    index0 = i*N*N + j*N + k;
                    xx[index] = x[index0];
                    yy[index] = y[index0];
                }

    fft_plan     = rfftw3d_create_plan( N, N, N, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE );
    fft_plan_inv = rfftw3d_create_plan( N, N, N, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE );

    rfftwnd_one_real_to_complex( fft_plan, xx, kxx );
    rfftwnd_one_real_to_complex( fft_plan, yy, kyy );

    for ( i=0; i<N; i++ )
        for( j=0; j<N; j++ )
                for( k=0; k<N3/2; k++ ) {
                    index = i*N*N3/2 + j*N3/2 + k;
                    a = kxx[index].re;
                    b = kxx[index].im;
                    kxx[index].re = a * kyy[index].re - b * kyy[index].im;
                    kxx[index].im = a * kyy[index].im + b * kyy[index].re;
                }

    rfftwnd_one_complex_to_real( fft_plan_inv, kxx, xx );

    rfftwnd_destroy_plan( fft_plan );
    rfftwnd_destroy_plan( fft_plan_inv );

/*
    for( i=0; i<10; i++ ) {
        printf( "[%i] %lf, %lf\n", i, x[i], (double)xx[i] / ( N*N*N ) );
    }
    for( i=0; i<10; i++ ) {
        printf( "[%i] %lf, %lf\n", N-i-2, x[N-i-2], (double)xx[N-i-2] / ( N*N*N ) );
    }
*/

    for( i=0; i<N; i++ )
        for( j=0; j<N; j++ )
                for( k=0; k<N; k++ ) {
                    index = i*N*N3 + j*N3 + k;
                    index0 = i*N*N + j*N + k;
                    c[index0] = xx[index] / CUBE(N);
                }
    myfree( xx );
    myfree( yy );

}

void corr_Tdiff_dens() {

    double z[2], *T, *Dens, dx, *T0, *Dens0, Tmin, Tmax, Densmin, Densmax,
          // Tmm[4], Densmm[4],
           *Tdiff, *Densdiff, Tdiffbar, Densbar,
           *TDensCorr, *TDensCorr1d, dlogr, rmax, rmin, r, *data, *num;
    int N, N2, N3, i, j, k, index, rN;
    MPI_Status status;
    long p;
    FILE *fd;

    writelog( "Correlation of T difference and density ...\n" );

    MPI_Gather( &All.RedShift, 1, MPI_DOUBLE, z, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if ( ThisTask == 0 ) {
        if ( fabs( z[0] - z[1] ) > 1e-5 ) {
            printf( "The same redshift is required by corralation of "
                "T difference and density\n" );
            endrun( 20190720 );
        }
    }

    N = All.NGrid;
    N2 = SQR(N);
    N3 = N2*N;
    dx = All.BoxSize / N;
    rN = All.CorrRN;
    rmin = dx * 2;
    rmax = All.BoxSize;
    dlogr = log( rmax/rmin ) / rN;

    mymalloc1( T, sizeof( double ) * N3 );
    mymalloc1( Dens, sizeof( double ) * N3 );
    mymalloc1( num, sizeof( double ) * N3 );
    mymalloc1( data, sizeof( double ) * N_Gas * 4 );

    if ( ThisTask == 0 ) {
        mymalloc1( T0, sizeof( double ) * N3 );
        mymalloc1( Dens0, sizeof( double ) * N3 );

        mymalloc2( Densdiff, sizeof( double ) * N3 );
        mymalloc2( Tdiff, sizeof( double ) * N3 );
        mymalloc2( TDensCorr, sizeof( double ) * N3 );
        mymalloc2( TDensCorr1d, sizeof( double ) * rN );
    }

    Tmin = Densmin = DBL_MAX;
    Tmax= Densmax = -DBL_MAX;

    for ( p=0; p<N_Gas; p++ ) {
        Tmin = ( SphP[p].Temp < Tmin ) ? SphP[p].Temp : Tmin;
        Tmax = ( SphP[p].Temp > Tmax ) ? SphP[p].Temp : Tmax;

        for( i=0; i<3; i++ )
            data[p+i] = P[p].Pos[i];
        data[p+3] = SphP[p].Temp;
    }

    field_to_grid( data, T, num, N_Gas, 1 );

    for ( p=0; p<N_Gas; p++ ) {
        Densmin = ( SphP[p].Density < Densmin ) ? SphP[p].Density : Densmin;
        Densmax = ( SphP[p].Density < Densmax ) ? SphP[p].Density : Densmax;
        for( i=0; i<3; i++ )
            data[p+i] = P[p].Pos[i];
        data[p+3] = SphP[p].Density;
    }

    field_to_grid( data, Dens, num, N_Gas, 1 );

    /*
    MPI_Gather( &Tmin, 1, MPI_DOUBLE, Tmm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Gather( &Tmax, 1, MPI_DOUBLE, Tmm+2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    t = Tmm[1];
    Tmm[1] = Tmm[2];
    Tmm[2] = t;

    MPI_Gather( &Densmin, 1, MPI_DOUBLE, Densmm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Gather( &Densmax, 1, MPI_DOUBLE, Densmm+2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    t = Densmm[1];
    Densmm[1] = Densmm[2];
    Densmm[2] = t;
    */


    if ( ThisTask == 0 ) {
        MPI_Recv( T0, N3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status  );
        MPI_Recv( Dens0, N3, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status  );
    }
    else {
        MPI_Send( T, N3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
        MPI_Send( Dens, N3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
    }
    //MPI_Sendrecv( T, N3, MPI_DOUBLE, 0, 0, T0, N3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status );
    //MPI_Sendrecv( Dens, N3, MPI_DOUBLE, 0, 11, Dens0, N3, MPI_DOUBLE, 1, 11, MPI_COMM_WORLD, &status );

    if ( ThisTask == 0 ) {

            Tdiffbar = Densbar = 0;
            for( index=0; index<N3; index++ ) {

                if ( T0[index] > 0 ) {
                    Tdiff[index] = ( T[index] - T0[index] ) / T0[index];
                }

                if ( Dens0[index] > 0 ) {
                    Densdiff[index] = ( Dens[index]-Dens0[index] ) / Dens0[index];
                }

                Tdiffbar += Tdiff[index];
                Densbar += Dens[index];

            }

            Tdiffbar /= N3;
            Densbar /= N3;

            FILE  *fd;
            fd = fopen( "T_Dens.dat", "w" );
            for( i=0; i<N; i++ )
                for( j=0; j<N; j++ )
                    for( k=0; k<N; k++ ) {
                        index = i * N2 + j * N + k;
                        fprintf( fd, "%g %g %g %g %g %g %g %g %g\n",
                            i * dx, j * dx, k * dx,
                            T[index], T0[index], Tdiff[index],
                            log10(Dens[index]), log10(Dens0[index]),
                            Densdiff[index]
                        );
            }
            fclose( fd );

            for( index=0; index<N3; index++ ) {
                Tdiff[index] = ( Tdiff[index] - Tdiffbar ) / Tdiffbar;
                Dens[index] = ( Dens[index] - Densbar ) / Densbar;
            }

        }

    if ( ThisTask == 0 ) {
        corr( Tdiff, Dens, TDensCorr, N );

        memset( num, 0, sizeof(int) * rN );
        for( i=0; i<N; i++ )
            for ( j=0; j<N; j++ )
                for ( k=0; k<N; k++ ) {
                    r = sqrt( i*i + j*j + k*k ) * dx;
                    p = log( r / rmin ) / dlogr; 
                    p = ( p < 0 ) ? 0 : p;
                    p = ( p > rN-1 ) ? rN-1 : p;
                    TDensCorr1d[p] += TDensCorr[ i * N2 + j * N + k ];
                    num[p] ++;
                }

        for( i=0; i<rN; i++ )
            if ( num[i] > 0 )
                TDensCorr1d[i] /= num[i];


        fd = fopen( "TDensCorr1d.dat", "w" );
        for( i=0; i<rN; i++ )
            fprintf( fd, "%g %e\n", rmin * exp( i*dlogr ), TDensCorr1d[i] );
        fclose( fd );

    }

    myfree( T );
    myfree( num );
    myfree( Dens );
    myfree( data );

    if ( ThisTask == 0 ) {
        myfree( T0 );
        myfree( Dens0 );
        myfree( Densdiff );
        myfree( Tdiff );
        myfree( TDensCorr );
        myfree( TDensCorr1d );
    }


}
