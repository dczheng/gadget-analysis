#include "allvars.h"

void c_analysis() {
    double *crc, crc_max, crc_min, log_crc_max, log_crc_min,
           dx, dy, x, y, glob_log_crc_max, glob_log_crc_min;
    int i, j, xi, yi;
    char buf[100];
    fprintf( LogFilefd, "cosmic ray c analysis ...\n" );
    crc = malloc( sizeof( double ) * para.PicSize * para.PicSize );
    memset( crc, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
    crc_max = -DBL_MAX;
    crc_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ){
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        crc[ xi*para.PicSize + yi ] += SphP[i].CR_C0;
    }

    for ( i=0; i<para.PicSize*para.PicSize; i++ ) {
        if ( crc[i] > 0 ) {
            if ( crc[i] > crc_max )
                crc_max = crc[i];
            if ( crc[i] < crc_min )
                crc_min = crc[i];
        }
    }
    if ( crc_max == -DBL_MAX )
        log_crc_max = -DBL_MAX;
    else
        log_crc_max = log10( crc_max );

    if ( crc_min == DBL_MAX )
        log_crc_min = DBL_MAX;
    else
        log_crc_min = log10( crc_min );

    fprintf( LogFilefd, "crc_max: %g, crc_min: %g\n"
            "log_crc_max: %g, log_crc_min: %g\n",
            crc_max, crc_min,
            log_crc_max, log_crc_min );

    if ( ThisTask == 0 )
    if ( access( "./crc/", 0 ) == -1 ){
        printf( "create directory ./crc by task %d\n", ThisTask);
        if ( mkdir( "./crc", 0755) == -1 ){
            printf( "failed create directory ./crc.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_crc_max, &glob_log_crc_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_crc_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_crc_min, &glob_log_crc_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_crc_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    fprintf( LogFilefd, "glob_log_crc_max: %g, glob_log_crc_min: %g\n",
            glob_log_crc_max, glob_log_crc_min );
    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( crc[i] > 0 )
            crc[i] = log10( crc[i] );
        else
            crc[i] = glob_log_crc_min;
    }
    sprintf( cb_label, "(10^x)" );
    sprintf( buf, "./crc/crc_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, crc, 0, para.PicSize, 0, para.PicSize,
            glob_log_crc_min, glob_log_crc_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "cosmic ray c (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_crc_min, glob_log_crc_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    fprintf( LogFilefd, "cosmic ray c analysis ... done.\n" );
}


void cosmic_ray_analysis() {
    c_analysis();
}
