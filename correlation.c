#include "allvars.h"
#include "drfftw.h"

double dr, rmax, rmin, r, dx, dx3;
int rN, NGrid, NGrid2;
long NGrid3;

void set_global_vars() {

    NGrid = All.NGrid;
    NGrid2 = SQR(NGrid);
    NGrid3 = NGrid2*NGrid;

    dx = BoxSize / NGrid;
    dx3 = CUBE( dx );
    rmin = dx * 2;
    rmax = BoxSize / 2.0;
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
    writelog( "compute correlation function ... done.\n" );

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

void output_grid_slice( double *grid, int N, char *fn_prefix ){
    // for test

    FILE *fd;
    int i, j, k, m;
    long index;

    for( m=0; m<3; m++ ) {
        fd = myfopen( "w", "%s_%03i_%c.dat", fn_prefix, SnapIndex, 'x'+m );
        for ( i=0; i<NGrid; i++ ) {
            for( j=0; j<NGrid; j++ ) {
                r = 0;
                for( k=0; k<N; k++ ) {
                    if ( m == 0 )
                        index =  k*NGrid2 + i*NGrid + j;
                    if ( m == 1 )
                        index =  j*NGrid2 + k*NGrid + i;
                    if ( m == 2 )
                        index =  i*NGrid2 + j*NGrid + k;
                    r += grid[index];
                }
                fprintf( fd, "%g ", r );
            }
            fprintf( fd, "\n" );
        }
        fclose( fd );
    }

}

void put_dm_on_grid( double *grid ) {

    long p;
    double *data;

    mymalloc1( data, sizeof( double ) * NumPart6[1] );

    for( p=0; p<NumPart6[1]; p++ )
        data[p] = 1;

    field_to_grid( data, grid, 1, 0, 0 );
    myfree( data );

}

void put_gas_on_grid( double *grid ) {

    long p;
    double *data;

    mymalloc1( data, sizeof( double ) * N_Gas );
    for( p=0; p<N_Gas; p++ )
        data[p] = P[p].Mass / dx3;

    field_to_grid( data, grid, 0, 0, 0 );

    output_grid_slice( grid, 10, "Dens" );
    myfree( data );

}

void put_dens_on_grid( double *grid, int Nmin ) {

    long p;
    double *data;

    mymalloc1( data, sizeof( double ) * N_Gas );
    for( p=0; p<N_Gas; p++ )
        data[p] = SphP[p].Density;

    field_to_grid( data, grid, 0, Nmin, 1 );

    output_grid_slice( grid, 10, "Dens" );
    myfree( data );

}

#ifdef PDFTDIFFDENS
void put_temp_on_grid( double *grid, int Nmin ) {

    long p;
    double *data;

    mymalloc1( data, sizeof( double ) * N_Gas );
    for( p=0; p<N_Gas; p++ )
        data[p] = SphP[p].Temp;

    field_to_grid( data, grid, 0, Nmin, 1 );
    myfree( data );

}
#endif

#ifdef CORRDM
void corr_dm() {

    double *DM, *DMCorr, *DMCorr1d, DMBar;
    long p, num_dm;
    int i;
    FILE *fd;

    writelog( "Correlation function of dark matter ...\n" );
    set_global_vars();

    num_dm = NumPart6[1];

    //printf( "%li, %li\n", offset_dm, num_dm );

    mymalloc1( DM, sizeof( double ) * NGrid3 );
    mymalloc1( DMCorr, sizeof( double ) * NGrid3 );
    mymalloc1( DMCorr1d, sizeof( double ) * rN );

    put_dm_on_grid( DM );

    DMBar = ((double)num_dm) / NGrid3;
    for( p=0; p<NGrid3; p++ ) {
        DM[p] = DM[p] / DMBar - 1;
    }

    corr( DM, DM, DMCorr );
    get_corr1d( DMCorr, DMCorr1d );

    fd = myfopen( "w", "DMCorr1d_%03i.dat", SnapIndex );
    for( i=0; i<rN; i++ )
            fprintf( fd, "%g %e\n", rmin +  i*dr, DMCorr1d[i] );
    fclose( fd );

    myfree( DM );
    myfree( DMCorr );
    myfree( DMCorr1d );

}
#endif

#ifdef CORRGAS
void corr_gas() {

    double *Gas, Gasbar, *GasCorr, *GasCorr1d;
    long p;
    int i;
    FILE *fd;

    writelog( "Correlation function of gas density ...\n" );
    set_global_vars();

    mymalloc1( Gas, sizeof( double ) * NGrid3 );
    mymalloc1( GasCorr, sizeof( double ) * NGrid3 );
    mymalloc1( GasCorr1d, sizeof( double ) * rN );

    put_gas_on_grid( Gas );

    Gasbar = 0;
    for( p=0; p<NGrid3; p++ ) {
        Gasbar += Gas[p];
    }
    Gasbar /= NGrid3;
    //printf( "%g %g\n", Gasbar, RhoBaryon );

    for ( p=0; p<NGrid3; p++ )
        Gas[p] = ( Gas[p] - Gasbar ) / Gasbar;

    corr( Gas, Gas, GasCorr );
    get_corr1d( GasCorr, GasCorr1d );

    fd = myfopen( "w", "GasCorr1d_%03i.dat", SnapIndex );
    for( i=0; i<rN; i++ )
        fprintf( fd, "%g %e\n", rmin +  i*dr, GasCorr1d[i] );
    fclose( fd );

    myfree( Gas );
    myfree( GasCorr );
    myfree( GasCorr1d );

}
#endif

#ifdef CORRTDIFFDENS
void corr_Tdiff_dens() {

}
#endif

#ifdef PDFTDIFFDENS
void pdf_Tdiff_dens() {

    double *Tdiff, *Dens, *Dens0, *tmp,
        Tdiffmax, Tdiffmin, Tdiffmin_1, Tdiffmax_1,
        dlogTdiff, dlogTdiff_1, dTdiff, dTdiff_1,
        Densmax, Densmin, *T, *T0, *imgT, *imgDens,
        Tmin, Tmax, dlogT,
        dlogDens, z[2];
    long p, Nmin;
    int N, Nhalf, N_p, N_m, i, j;
    MPI_Status status;

    MPI_Gather( &Redshift, 1, MPI_DOUBLE, z, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if ( ThisTask == 0 )
        if ( fabs( z[0] - z[1] ) > 1e-5 ) {
            printf( "The same redshift is required by corralation of "
                "T difference and density\n" );
            endrun( 20190720 );
        }

    writelog( "PDF of Tdiff and denssity ... \n" );
    set_global_vars();

    mymalloc1( Dens, sizeof(double)*NGrid3 );
    mymalloc1( T, sizeof(double)*NGrid3 );

    if ( ThisTask == 0 ) {
        mymalloc1( Tdiff, sizeof(double)*NGrid3 );
        mymalloc1( Dens0, sizeof(double)*NGrid3 );
        mymalloc1( T0, sizeof(double)*NGrid3 );
    }

    N = All.PicSize;   
    Nhalf = N / 2;
    N_p = N - Nhalf - 1;
    N_m = Nhalf;
    Nmin = 30;

    put_dens_on_grid( Dens, Nmin );
    put_temp_on_grid( T, Nmin );

    for( p=0; p<NGrid3; p++ )
        Dens[p] /= RhoBaryon;

    if ( ThisTask == 0 ) {
        MPI_Recv( Dens0, NGrid3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status  );
        MPI_Recv( T0, NGrid3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status  );
    }
    else{
        MPI_Send( Dens, NGrid3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
        MPI_Send( T, NGrid3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
    }

    if ( ThisTask == 0 ) {

        for( i=0; i<NGrid3; i++ ) {
            if ( T[i] == 0 || T0[i] == 0 )
                Tdiff[i] = 0;
            else
                Tdiff[i] = (T[i] - T0[i]) / T0[i];
            //Tdiff[i] = T[i] - T0[i];

            if ( fabs(Tdiff[i]) < 1e-3 )
                Tdiff[i] = 0;
        }

        //tmp = Dens;
        //Dens = T;
        //T = tmp;

        Tmax = Tdiffmax_1 = Tdiffmax = Densmax = -DBL_MAX;
        Tmin = Tdiffmin_1 = Tdiffmin = Densmin = DBL_MAX;

        for( p=0; p<NGrid3; p++ ) {

            vmax2( Densmax, Dens[p] );
            vmin20( Densmin, Dens[p] );

            vmax2( Tmax, T[p] );
            vmin20( Tmin, T[p] );


            if ( Tdiff[p] > 0 ) {
                vmax2( Tdiffmax, Tdiff[p] );
                vmin2( Tdiffmin, Tdiff[p] );
            }

            if ( Tdiff[p] < 0 ) {
                vmax2( Tdiffmax_1, Tdiff[p] );
                vmin2( Tdiffmin_1, Tdiff[p] );
            }
        }
    }

    // Those sync is only for paraller output
    MPI_Bcast( &Tdiffmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &Tdiffmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &Tdiffmin_1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &Tdiffmax_1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &Densmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &Densmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &Tmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &Tmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    writelog( "Tdiff, min: %g, max: %g, min_1: %g, max_1: %g\n",
          Tdiffmin, Tdiffmax, Tdiffmin_1, Tdiffmax_1 );
    writelog( "Dens,  min: %g, max: %g\n", Densmin,  Densmax );
    writelog( "T,  min: %g, max: %g\n", Tmin,  Tmax );

    //Tdiffmin   =    0.001;
    //Tdiffmax_1 = -0.001;


    if ( ThisTask == 0 ) {

        /*
        Tdiffmax   = 1;
        Tdiffmin_1 = -1;

        double n1, n2;
        n1 = n2 = 0;
        for ( p=0; p<NGrid3; p++ ) {
            if ( Tdiff[p] > 1 )
                n1++;
            if ( Tdiff[p] < -1 )
                n2++;
        }
        printf( "n1: %i, n2: %i\n", n1, n2 );
        */

        //Tdiffmax = 10;
        //Tdiffmin_1 = -10;

        dlogDens = log10( Densmax / Densmin ) / N;
        dlogT = log10( Tmax / Tmin ) / N;

        dlogTdiff = log10( Tdiffmax / Tdiffmin ) / N_p; 
        dlogTdiff_1 = log10( Tdiffmax_1 / Tdiffmin_1 ) / N_m; 

        dTdiff = ( Tdiffmax - Tdiffmin ) / N_p; 
        dTdiff_1 = ( Tdiffmax_1 - Tdiffmin_1 ) / N_m; 

        reset_img();

        mymalloc2( imgT, sizeof(double) * N * N );
        mymalloc2( imgDens, sizeof(double) * N * N );

        for( p=0; p<NGrid3; p++ ) {

            if ( Tdiff[p] < Tdiffmin && Tdiff[p] > Tdiffmax_1 )
                //i = Nhalf;
                continue;
            else {

                if ( Tdiff[p] > 0 )
                    i = log10( Tdiff[p] / Tdiffmin ) / dlogTdiff + Nhalf + 1;
                    //i = ( Tdiff[p] - Tdiffmin ) / dTdiff + Nhalf + 1;

                if ( Tdiff[p] < 0 )
                    i = log10( Tdiff[p] / Tdiffmin_1 ) / dlogTdiff_1;
                    //i = ( Tdiff[p] - Tdiffmin_1 ) / dTdiff_1;

            }

            /*
            if ( i == Nhalf ) {
                printf( "%g\n", Tdiff[p] );
                endrun( 20190801 );
            }
            */

            if ( i > N-1 || i < 0 ) {
                continue;
            }

            if ( Dens[p] == 0  )
                continue;
            j = log10( Dens[p] / Densmin ) / dlogDens;
            imgDens[ i*N + j ] ++;

            if ( T[p] == 0  )
                continue;
            j = log10( T[p] / Tmin ) / dlogT;
            imgT[ i*N + j ] ++;

        }

/*
        for ( i=0; i<N*N; i++ ) {
            if ( image.img[i] < 5 )
                image.img[i] = 0;
        }
*/

        img_xmin = Densmin; 
        img_xmax = Densmax; 

        img_ymin = Tdiffmin_1; 
        img_ymax = Tdiffmax; 

        img_props(0) = Tdiffmin_1;
        img_props(1) = Tdiffmax_1;
        img_props(2) = dlogTdiff_1;
        img_props(3) = Tdiffmin;
        img_props(4) = Tdiffmax;
        img_props(5) = dlogTdiff;
        img_props(6) = Densmin;
        img_props(7) = Densmax;
        img_props(8) = dlogDens;

        image.img =imgDens; 
        write_img( "PdfTdiffDens.dat", "PpdfTdiffDens" );

        img_xmin = Densmin; 
        img_xmax = Densmax; 

        img_props(6) = Tmin;
        img_props(7) = Tmax;
        img_props(8) = dlogT;

        image.img =imgT; 
        write_img( "PdfTdiffT.dat", "PpdfTdiffT" );

        myfree( imgT );
        myfree( imgDens );
        myfree( Tdiff );
        myfree( Dens0 );
        myfree( T0 );
    }

    myfree( Dens );
    myfree( T );

    do_sync( "" );

}
#endif
