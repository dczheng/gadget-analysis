#include "allvars.h"

void hg_electrons_analysis() {
    double *rho_n, x, y, dx, dy, rho_n_max, rho_n_min,
           log_rho_n_min, log_rho_n_max,
           glob_log_rho_n_min, glob_log_rho_n_max;
    char buf[100];
    int i,j, xi, yi, zi;
    fprintf( LogFilefd, "high energy electrons analysis...\n" );
    rho_n = ( double* ) malloc( sizeof(double) * para.PicSize * para.PicSize );
    memset( rho_n, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = BoxSize / para.PicSize;
    dy = BoxSize / para.PicSize;

    for ( i=0; i<N_Gas; i++ ) {
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        rho_n[xi*para.PicSize+yi] += SphP[i].CRE_n0 * SphP[i].Density * cgs_me / ( cgs_g/CUBE(cgs_cm) );
        //printf( "xi: %d, yi: %d, %g\n", xi, yi, rho_n[xi*para.PicSize+yi] );
    }
    rho_n_max = -DBL_MAX;
    rho_n_min = DBL_MAX;
    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( rho_n[i] > 0 ){
            if ( rho_n[i] < rho_n_min )
                rho_n_min = rho_n[i];
            if ( rho_n[i] > rho_n_max )
                rho_n_max = rho_n[i];
        }
    }
    log_rho_n_max = log10( rho_n_max );
    log_rho_n_min = log10( rho_n_min );
    fprintf( LogFilefd, "rho_n_max: %g, rho_n_min: %g \n"
             "log_rho_n_max: %g log_rho_n_min: %g\n",
             rho_n_max, rho_n_min, log_rho_n_max, log_rho_n_min );

    if ( ThisTask == 0 )
    if ( access( "./hge_n/", 0 ) == -1 ){
        printf( "create directory ./hge_n/ by task %d\n", ThisTask );
        if ( mkdir( "./hge_n/", 0755) == -1 ){
            printf( "failed create directory ./hge_n/.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_rho_n_max, &glob_log_rho_n_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_rho_n_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_rho_n_min, &glob_log_rho_n_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_rho_n_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    fprintf( LogFilefd, "glob_log_rho_n_max: %g, glob_log_rho_n_min: %g\n",
            glob_log_rho_n_max, glob_log_rho_n_min );
    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( rho_n[i] > 0 )
            rho_n[i] = log10( rho_n[i] );
        else
            rho_n[i] = glob_log_rho_n_min;
    }
    sprintf( cb_label, "(10^x)" );
    sprintf( buf, "./hge_n/hge_n_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, rho_n, 0, para.PicSize, 0, para.PicSize,
            glob_log_rho_n_min, glob_log_rho_n_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "high energy electron number densit (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_rho_n_min, glob_log_rho_n_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    free( rho_n );
    fprintf( LogFilefd, "high energy electrons analysis...done.\n" );
}

