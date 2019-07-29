#include "allvars.h"
#include "drfftw.h"

double dr, rmax, rmin, r, *data, dx;
int rN, NGrid, NGrid2;
long NGrid3;

void set_global_vars() {

    NGrid = All.NGrid;
    NGrid2 = SQR(NGrid);
    NGrid3 = NGrid2*NGrid;

    dx = All.BoxSize / NGrid;
    rmin = dx * 2;
    rmax = All.BoxSize / 2.0;
    dr = dx * 2;
    rN = ( rmax - rmin ) / dr;

}

void corr( double *x, double *y, double *c ) {

    /*
    \int x(t) + y(t+r) dt  = \int x(-t) + y(r-t) dt
    */

    long i, j, k, index, index0, FFT_NGrid3;
    double a, b, fac;

    fftw_real *xx, *yy;
    fftw_complex *kxx, *kyy;
    rfftwnd_plan fft_plan, fft_plan_inv;
    FFT_NGrid3 = 2 * ( NGrid/2+1 );

    writelog( "compute correlation function ...\n" );
    mymalloc2( xx, NGrid2*FFT_NGrid3*sizeof( fftw_real ) );
    mymalloc2( yy, NGrid2*FFT_NGrid3*sizeof( fftw_real ) );
    kxx = ( fftw_complex* ) xx;
    kyy = ( fftw_complex* ) yy;
    fac = SQR(1.0/NGrid3);

    //printf( "%i %i %li %li\n", NGrid, NGrid2, NGrid3, FFT_NGrid3 );
    for( i=0; i<NGrid; i++ )
        for( j=0; j<NGrid; j++ )
                for( k=0; k<NGrid; k++ ) {
                    index = i*NGrid*FFT_NGrid3 + j*FFT_NGrid3 + k;

                    index0 = (NGrid-1-i)*NGrid2 + (NGrid-1-j)*NGrid + (NGrid-1-k);
                    xx[index] = x[index0];

                    index0 = i*NGrid2 + j*NGrid + k;
                    yy[index] = y[index0];
                }

    fft_plan     = rfftw3d_create_plan( NGrid, NGrid, NGrid, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE );
    fft_plan_inv = rfftw3d_create_plan( NGrid, NGrid, NGrid, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE );

    rfftwnd_one_real_to_complex( fft_plan, xx, NULL );
    rfftwnd_one_real_to_complex( fft_plan, yy, NULL );

    for ( i=0; i<NGrid; i++ )
        for( j=0; j<NGrid; j++ )
                for( k=0; k<FFT_NGrid3/2; k++ ) {
                    index = i*NGrid*FFT_NGrid3/2 + j*FFT_NGrid3/2 + k;
                    a = kxx[index].re;
                    b = kxx[index].im;
                    kxx[index].re = a * kyy[index].re - b * kyy[index].im;
                    kxx[index].im = a * kyy[index].im + b * kyy[index].re;
                }

    rfftwnd_one_complex_to_real( fft_plan_inv, kxx, NULL );

    rfftwnd_destroy_plan( fft_plan );
    rfftwnd_destroy_plan( fft_plan_inv );

    for( i=0; i<NGrid; i++ )
        for( j=0; j<NGrid; j++ )
                for( k=0; k<NGrid; k++ ) {
                    index = i*NGrid*FFT_NGrid3 + j*FFT_NGrid3 + k;
                    index0 = i*NGrid2 + j*NGrid + k;
                    c[index0] = xx[index] * fac;
                }
    myfree( xx );
    myfree( yy );
    writelog( "compute correlation function... done.\n" );

}

void get_corr1d( double *corr3d, double *corr1d ) {

    int i, j, k, ii, jj, kk, index, *num;
    
    mymalloc2( num, sizeof(int) * rN );
    for( i=0; i<NGrid; i++ )
        for ( j=0; j<NGrid; j++ )
            for ( k=0; k<NGrid; k++ ) {

                ii = ( i > NGrid/2 ) ? i-NGrid : i;
                jj = ( j > NGrid/2 ) ? j-NGrid : j;
                kk = ( k > NGrid/2 ) ? k-NGrid : k;

                r = sqrt( ii*ii + jj*jj + kk*kk ) * dx;

                index = (r-rmin) / dr + 1; 

                if ( index<0 || index > rN-1 )
                    continue;

                corr1d[index] += corr3d[ i * NGrid2 + j * NGrid + k ];
                num[index] ++;
            }

    for( i=0; i<rN; i++ )
        if ( num[i] > 0 )
           corr1d[i] /= num[i];

    myfree( num );

}

void corr_dm() {

    double *DM, *DMCorr, *DMCorr1d, DMBar;
    long p, num_dm, offset_dm;
    int i;
    FILE *fd;
    char buf[100];

    writelog( "Correlation function of dark matter ...\n" );
    set_global_vars();

    offset_dm = OffsetPart6[1];
    num_dm = NumPart6[1];

    //printf( "%li, %li\n", offset_dm, num_dm );

    mymalloc1( data, sizeof( double ) * num_dm * 4 );
    mymalloc1( DM, sizeof( double ) * NGrid3 );

    mymalloc1( DMCorr, sizeof( double ) * NGrid3 );
    mymalloc1( DMCorr1d, sizeof( double ) * rN );

    for( p=0; p<num_dm; p++ ) {
        for( i=0; i<3; i++ )
            data[4*p+i] = P[offset_dm+p].Pos[i];
        data[4*p+3] = 1;
    }

    field_to_grid( data, DM, num_dm, 0 );

/*
    fd = fopen( "dm.dat", "w" );
    for ( i=0; i<NGrid; i++ ) {
        for( j=0; j<NGrid; j++ ) {
            r = 0;
            for( k=0; k<100; k++ )
                r += DM[ i * NGrid2 + j * NGrid + k ];
            fprintf( fd, "%g ", r );
        }
        fprintf( fd, "\n" );
    }
    fclose( fd );
*/

    DMBar = ((double)num_dm) / NGrid3;
    for( p=0; p<NGrid3; p++ ) {
        DM[p] = DM[p] / DMBar - 1;
    }

    corr( DM, DM, DMCorr );
    get_corr1d( DMCorr, DMCorr1d );

/*
    fd = fopen( "dm_corr3d.dat", "w" );
    for ( i=0; i<NGrid; i++ ) {
        for( j=0; j<NGrid; j++ ) {
            r = 0;
            for( k=0; k<10; k++ )
                r += DMCorr[ i * NGrid2 + j * NGrid + k ];
            fprintf( fd, "%g ", r );
        }
        fprintf( fd, "\n" );
    }
    fclose( fd );
*/

    sprintf( buf, "DMCorr1d_%03i.dat", All.SnapIndex );
    fd = fopen( buf, "w" );
    for( i=0; i<rN; i++ )
            fprintf( fd, "%g %e\n", rmin +  i*dr, DMCorr1d[i] );
    fclose( fd );

    myfree( data );
    myfree( DM );
    myfree( DMCorr );
    myfree( DMCorr1d );
    put_sep;

}

void corr_dens() {

    double *Dens, Densbar, *DensCorr, *DensCorr1d;
    long p;
    int i;
    FILE *fd;
    char buf[100];

    writelog( "Correlation function of gas density ...\n" );
    set_global_vars();

    mymalloc1( data, sizeof( double ) * N_Gas * 4 );
    mymalloc1( Dens, sizeof( double ) * NGrid3 );
    mymalloc1( DensCorr, sizeof( double ) * NGrid3 );
    mymalloc1( DensCorr1d, sizeof( double ) * rN );

    Densbar = 0;
    for( p=0; p<N_Gas; p++ ) {
        for( i=0; i<3; i++ )
            data[4*p+i] = P[p].Pos[i];
        data[4*p+3] = SphP[p].Density;
        Densbar += SphP[p].Density;
    }
    Densbar /= N_Gas;
    field_to_grid( data, Dens, N_Gas, 1 );

    for ( p=0; p<NGrid3; p++ )
        Dens[p] = ( Dens[p] - Densbar ) / Densbar;

    corr( Dens, Dens, DensCorr );
    get_corr1d( DensCorr, DensCorr1d );

    sprintf( buf, "DensCorr1d_%03i.dat", All.SnapIndex );
    fd = fopen( buf, "w" );
    for( i=0; i<rN; i++ )
        fprintf( fd, "%g %e\n", rmin +  i*dr, DensCorr1d[i] );
    fclose( fd );

    myfree( data );
    myfree( Dens );
    myfree( DensCorr );
    myfree( DensCorr1d );
    put_sep;

}

void corr_Tdiff_dens() {

    double z[2], *T, *Dens, *T0, *Dens0, Tmin, Tmax, Densmin, Densmax,
          // Tmm[4], Densmm[4],
           *Tdiff, *Densdiff, Tdiffbar, Densbar,
           *TDensCorr, *TDensCorr1d;
    int i, index;
    MPI_Status status;
    long p;
    FILE *fd;

    writelog( "Correlation function of T difference and density ...\n" );
    set_global_vars();

    MPI_Gather( &All.RedShift, 1, MPI_DOUBLE, z, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if ( ThisTask == 0 ) {
        if ( fabs( z[0] - z[1] ) > 1e-5 ) {
            printf( "The same redshift is required by corralation of "
                "T difference and density\n" );
            endrun( 20190720 );
        }
    }

    mymalloc1( T, sizeof( double ) * NGrid3 );
    mymalloc1( Dens, sizeof( double ) * NGrid3 );
    mymalloc1( data, sizeof( double ) * N_Gas * 4 );

    if ( ThisTask == 0 ) {
        mymalloc1( T0, sizeof( double ) * NGrid3 );
        mymalloc1( Dens0, sizeof( double ) * NGrid3 );

        mymalloc2( Densdiff, sizeof( double ) * NGrid3 );
        mymalloc2( Tdiff, sizeof( double ) * NGrid3 );
        mymalloc2( TDensCorr, sizeof( double ) * NGrid3 );
        mymalloc2( TDensCorr1d, sizeof( double ) * rN );
    }

    Tmin = Densmin = DBL_MAX;
    Tmax= Densmax = -DBL_MAX;

    for( p=0; p<N_Gas; p++ ) {
        for( i=0; i<3; i++ )
            data[4*p+i] = P[p].Pos[i];
        data[4*p+3] = SphP[p].Temp;

        Tmin = ( SphP[p].Temp < Tmin ) ? SphP[p].Temp : Tmin;
        Tmax = ( SphP[p].Temp > Tmax ) ? SphP[p].Temp : Tmax;

    }

    field_to_grid( data, T, N_Gas, 1 );

    for( p=0; p<N_Gas; p++ ) {
        for( i=0; i<3; i++ )
            data[4*p+i] = P[p].Pos[i];
        data[4*p+3] = SphP[p].Density;

        Densmin = ( SphP[p].Density < Densmin ) ? SphP[p].Density : Densmin;
        Densmax = ( SphP[p].Density < Densmax ) ? SphP[p].Density : Densmax;

    }

    field_to_grid( data, Dens, N_Gas, 1 );

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
        MPI_Recv( T0, NGrid3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status  );
        MPI_Recv( Dens0, NGrid3, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status  );
    }
    else {
        MPI_Send( T, NGrid3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
        MPI_Send( Dens, NGrid3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD );
    }
    //MPI_Sendrecv( T, NGrid3, MPI_DOUBLE, 0, 0, T0, NGrid3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status );
    //MPI_Sendrecv( Dens, NGrid3, MPI_DOUBLE, 0, 11, Dens0, NGrid3, MPI_DOUBLE, 1, 11, MPI_COMM_WORLD, &status );

    if ( ThisTask == 0 ) {

            Tdiffbar = Densbar = 0;
            for( index=0; index<NGrid3; index++ ) {

                if ( T0[index] > 0 ) {
                    Tdiff[index] = ( T[index] - T0[index] ) / T0[index];
                }

                if ( Dens0[index] > 0 ) {
                    Densdiff[index] = ( Dens[index]-Dens0[index] ) / Dens0[index];
                }

                Tdiffbar += Tdiff[index];
                Densbar += Dens[index];

            }

            Tdiffbar /= NGrid3;
            Densbar /= NGrid3;

            /*
            FILE  *fd;
            fd = fopen( "T_Dens.dat", "w" );
            for( i=0; i<N; i++ )
                for( j=0; j<N; j++ )
                    for( k=0; k<N; k++ ) {
                        index = i * NGrid2 + j * N + k;
                        fprintf( fd, "%g %g %g %g %g %g %g %g %g\n",
                            i * dx, j * dx, k * dx,
                            T[index], T0[index], Tdiff[index],
                            log10(Dens[index]), log10(Dens0[index]),
                            Densdiff[index]
                        );
            }
            fclose( fd );
            */

            /*
            fd = fopen( "Dens.dat", "w" );
            for( p=0; p<N_Gas; p++ ) {
                fprintf( fd, "%g %g %g %g\n",
                        P[p].Pos[0],
                        P[p].Pos[1],
                        P[p].Pos[2],
                        log10(SphP[p].Density)
                        );
            }
            fclose( fd );
            */

        for( index=0; index<NGrid3; index++ ) {
            Tdiff[index] = ( Tdiff[index] - Tdiffbar ) / Tdiffbar;
            Dens[index] = ( Dens[index] - Densbar ) / Densbar;
        }

        corr( Tdiff, Dens, TDensCorr);
        get_corr1d( TDensCorr, TDensCorr1d );

        fd = fopen( "TDensCorr1d.dat", "w" );
        for( i=0; i<rN; i++ )
            fprintf( fd, "%g %e\n", rmin +  i*dr, TDensCorr1d[i] );
        fclose( fd );

    }

    myfree( T );
    myfree( Dens );

    if ( ThisTask == 0 ) {
        myfree( T0 );
        myfree( Dens0 );
        myfree( Densdiff );
        myfree( Tdiff );
        myfree( TDensCorr );
        myfree( TDensCorr1d );
    }

    put_sep;


}
