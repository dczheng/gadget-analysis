#include "allvars.h"

void mach_analysis(){
    int i,j, xi, yi, pt;
    double *mn, x, y, dx, dy, mn_max, mn_min, log_mn_min, log_mn_max,
           glob_log_mn_min, glob_log_mn_max;
    char buf[100];
    pt = 0;
    fprintf( LogFilefd, "mach number analysis ...\n" );
    mn = ( double* ) malloc( sizeof(double) * para.PicSize * para.PicSize );
    memset( mn, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
    mn_max = -DBL_MAX;
    mn_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        mn[ xi*para.PicSize + yi ] += SphP[i].MachNumber;
        if (SphP[i].MachNumber > mn_max)
            mn_max = SphP[i].MachNumber;
        if (SphP[i].MachNumber < mn_min)
            mn_min = SphP[i].MachNumber;
    }
    fprintf( LogFilefd, "SphP Max MachNumber: %g\n", mn_max );
    fprintf( LogFilefd, "SphP Min MachNumber: %g\n", mn_min );
    mn_max = -DBL_MAX;
    mn_min = DBL_MAX;
    for ( i=0; i<para.PicSize*para.PicSize; i++ ) {
        if ( mn[i] > 0 ){
            if ( mn[i] > mn_max )
                mn_max = mn[i];
            if ( mn[i] < mn_min )
                mn_min = mn[i];
        }
    }
    if ( mn_max == -DBL_MAX )
        log_mn_max = -DBL_MAX;
    else
    log_mn_max = log10( mn_max );
    if ( mn_min == DBL_MAX )
        log_mn_min = DBL_MAX;
    else
        log_mn_min = log10( mn_min );
    fprintf( LogFilefd, "mn_max: %g, mn_min: %g\n"
            "log_mn_max: %g, log_mn_min: %g\n",
            mn_max, mn_min, log_mn_max, log_mn_min );

    if ( ThisTask == 0 )
    if ( access( "./mn/", 0 ) == -1 ){
        printf( "create directory ./mn by task %d\n", ThisTask );
        if ( mkdir( "./mn/", 0755) == -1 ){
            printf( "failed create directory ./mn.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_mn_max, &glob_log_mn_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_mn_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_mn_min, &glob_log_mn_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_mn_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    for ( i=0; i<para.PicSize*para.PicSize; i++ )
        if ( mn[i] > 0 )
            mn[i] = log10( mn[i] );
        else
            mn[i] = glob_log_mn_min;
    fprintf( LogFilefd, "glob_log_mn_max: %g, glob_log_mn_min: %g\n",
            glob_log_mn_max, glob_log_mn_min );
    sprintf( cb_label, "(10^x)" );
    sprintf( buf, "./mn/mn_%.2f.png", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, mn, 0, para.PicSize, 0, para.PicSize,
            glob_log_mn_min, glob_log_mn_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "mach number (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_mn_min, glob_log_mn_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    free( mn );
    fprintf( LogFilefd, "mach number analysis ...done.\n" );
}

