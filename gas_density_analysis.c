#include "allvars.h"
void gas_density_analysis(){
    int i,j, xi, yi;
    double *rho, x, y, dx, dy, rho_max, rho_min;
    double log_rho_max, log_rho_min;
    double glob_log_rho_min, glob_log_rho_max;
    char buf[200];
    fprintf( LogFilefd, "gas density analysis ...\n");
    rho = ( double* ) malloc( sizeof(double) * para.PicSize * para.PicSize );
    memset( rho, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
    rho_max = -DBL_MAX;
    rho_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        rho[ xi*para.PicSize + yi ] += SphP[i].Density; // ( g/(cm*cm*cm) );
    }
    for ( i=0; i<para.PicSize*para.PicSize; i++ ) {
        if ( rho[i] > 0 ) {
            if ( rho[i] > rho_max )
                rho_max = rho[i];
            if ( rho[i] < rho_min )
                rho_min = rho[i];
        }
    }
    log_rho_max = log10( rho_max );
    log_rho_min = log10( rho_min );

    fprintf( LogFilefd, "rho_max: %g, rho_min: %g\n"
            "log_rho_max: %g, log_rho_min: %g\n",
            rho_max, rho_min, log_rho_max, log_rho_min );
    if ( ThisTask == 0 )
    if ( access( "./gas_rho/", 0 ) == -1 ){
        printf( "create directory ./gas_rho by task %d\n", ThisTask );
        if ( mkdir( "./gas_rho", 0755) == -1 ){
            printf( "failed create directory ./gas_rho.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_rho_max, &glob_log_rho_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_rho_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_rho_min, &glob_log_rho_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_rho_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    fprintf( LogFilefd, "glob_log_rho_max: %g, glob_log_rho_min: %g\n",
            glob_log_rho_max, glob_log_rho_min );
    for ( i=0; i<para.PicSize*para.PicSize; i++ )
        if ( rho[i] > 0 )
            rho[i] = log10( rho[i] );
        else
            rho[i] = glob_log_rho_min;

    sprintf( cb_label, "(10^x)   g/cm^3" );
    sprintf( buf, "./gas_rho/gas_rho_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, rho, 0, para.PicSize, 0, para.PicSize,
            glob_log_rho_min, glob_log_rho_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "gas density (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_rho_min, glob_log_rho_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    free( rho );
    fprintf( LogFilefd, "gas density analysis ...done.\n");
}

