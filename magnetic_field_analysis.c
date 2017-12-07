#include "allvars.h"

void B_analysis(){
    double *mag, mag_max, mag_min, log_mag_max, log_mag_min,
           dx, dy, x, y, glob_log_mag_max, glob_log_mag_min;
    int i, j, xi, yi;
    char buf[100];
    fprintf( LogFilefd, "magnetic field analysis ...\n" );
    mag = malloc( sizeof( double ) * para.PicSize * para.PicSize );
    memset( mag, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
    mag_max = -DBL_MAX;
    mag_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ){
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        mag[ xi*para.PicSize + yi ] += sqrt( pow( SphP[i].B[0], 2.0 ) +
                pow( SphP[i].B[1],2.0 ) + pow( SphP[i].B[2], 2.0 ) );
    }

    for ( i=0; i<para.PicSize*para.PicSize; i++ ) {
        if ( mag[i] > 0 ) {
            if ( mag[i] > mag_max )
                mag_max = mag[i];
            if ( mag[i] < mag_min )
                mag_min = mag[i];
        }
    }
    if ( mag_max == -DBL_MAX )
        log_mag_max = -DBL_MAX;
    else
        log_mag_max = log10( mag_max );

    if ( mag_min == DBL_MAX )
        log_mag_min = DBL_MAX;
    else
        log_mag_min = log10( mag_min );

    fprintf( LogFilefd, "mag_max: %g, mag_min: %g\n"
            "log_mag_max: %g, log_mag_min: %g\n",
            mag_max, mag_min,
            log_mag_max, log_mag_min );

    if ( ThisTask == 0 )
    if ( access( "./mag/", 0 ) == -1 ){
        printf( "create directory ./mag by task %d\n", ThisTask);
        if ( mkdir( "./mag", 0755) == -1 ){
            printf( "failed create directory ./mag.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_mag_max, &glob_log_mag_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_mag_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_mag_min, &glob_log_mag_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_mag_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    fprintf( LogFilefd, "glob_log_mag_max: %g, glob_log_mag_min: %g\n",
            glob_log_mag_max, glob_log_mag_min );
    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( mag[i] > 0 )
            mag[i] = log10( mag[i] );
        else
            mag[i] = glob_log_mag_min;
    }
    sprintf( cb_label, "(10^x)  Gauss" );
    sprintf( buf, "./mag/mag_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, mag, 0, para.PicSize, 0, para.PicSize,
            glob_log_mag_min, glob_log_mag_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "magnetic field (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_mag_min, glob_log_mag_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    fprintf( LogFilefd, "magnetic field analysis ... done.\n" );
}

void divB_analysis() {
    double *divB, divB_max, divB_min, log_divB_max, log_divB_min,
           dx, dy, x, y, glob_log_divB_max, glob_log_divB_min;
    int i, j, xi, yi;
    char buf[100];
    fprintf( LogFilefd, "divergence of magnetic field analysis ...\n" );
    divB = malloc( sizeof( double ) * para.PicSize * para.PicSize );
    memset( divB, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
    divB_max = -DBL_MAX;
    divB_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ){
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        divB[ xi*para.PicSize + yi ] += SphP[i].divB;
    }

    for ( i=0; i<para.PicSize*para.PicSize; i++ ) {
        if ( divB[i] > 0 ) {
            if ( divB[i] > divB_max )
                divB_max = divB[i];
            if ( divB[i] < divB_min )
                divB_min = divB[i];
        }
    }
    if ( divB_max == -DBL_MAX )
        log_divB_max = -DBL_MAX;
    else
        log_divB_max = log10( divB_max );

    if ( divB_min == DBL_MAX )
        log_divB_min = DBL_MAX;
    else
        log_divB_min = log10( divB_min );

    fprintf( LogFilefd, "divB_max: %g, divB_min: %g\n"
            "log_divB_max: %g, log_divB_min: %g\n",
            divB_max, divB_min,
            log_divB_max, log_divB_min );

    if ( ThisTask == 0 )
    if ( access( "./divB/", 0 ) == -1 ){
        printf( "create directory ./divB by task %d\n", ThisTask);
        if ( mkdir( "./divB", 0755) == -1 ){
            printf( "failed create directory ./divB.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_divB_max, &glob_log_divB_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_divB_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_divB_min, &glob_log_divB_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_divB_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    fprintf( LogFilefd, "glob_log_divB_max: %g, glob_log_divB_min: %g\n",
            glob_log_divB_max, glob_log_divB_min );
    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( divB[i] > 0 )
            divB[i] = log10( divB[i] );
        else
            divB[i] = glob_log_divB_min;
    }
    sprintf( cb_label, "(10^x)" );
    sprintf( buf, "./divB/divB_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, divB, 0, para.PicSize, 0, para.PicSize,
            glob_log_divB_min, glob_log_divB_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "divB (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_divB_min, glob_log_divB_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    fprintf( LogFilefd, "divergence of magnetic field analysis ... done.\n" );
}

void dBdt_analysis() {
    double *dBdt, dBdt_max, dBdt_min, log_dBdt_max, log_dBdt_min,
           dx, dy, x, y, glob_log_dBdt_max, glob_log_dBdt_min;
    int i, j, xi, yi;
    char buf[100];
    fprintf( LogFilefd, "rate of change magnetic field analysis ...\n" );
    dBdt = malloc( sizeof( double ) * para.PicSize * para.PicSize );
    memset( dBdt, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
    dBdt_max = -DBL_MAX;
    dBdt_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ){
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        dBdt[ xi*para.PicSize + yi ] += SphP[i].dBdt;
    }

    for ( i=0; i<para.PicSize*para.PicSize; i++ ) {
        if ( dBdt[i] > 0 ) {
            if ( dBdt[i] > dBdt_max )
                dBdt_max = dBdt[i];
            if ( dBdt[i] < dBdt_min )
                dBdt_min = dBdt[i];
        }
    }
    if ( dBdt_max == -DBL_MAX )
        log_dBdt_max = -DBL_MAX;
    else
        log_dBdt_max = log10( dBdt_max );

    if ( dBdt_min == DBL_MAX )
        log_dBdt_min = DBL_MAX;
    else
        log_dBdt_min = log10( dBdt_min );

    fprintf( LogFilefd, "dBdt_max: %g, dBdt_min: %g\n"
            "log_dBdt_max: %g, log_dBdt_min: %g\n",
            dBdt_max, dBdt_min,
            log_dBdt_max, log_dBdt_min );

    if ( ThisTask == 0 )
    if ( access( "./dBdt/", 0 ) == -1 ){
        printf( "create directory ./dBdt by task %d\n", ThisTask);
        if ( mkdir( "./dBdt", 0755) == -1 ){
            printf( "failed create directory ./dBdt.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_dBdt_max, &glob_log_dBdt_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_dBdt_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_dBdt_min, &glob_log_dBdt_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_dBdt_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    fprintf( LogFilefd, "glob_log_dBdt_max: %g, glob_log_dBdt_min: %g\n",
            glob_log_dBdt_max, glob_log_dBdt_min );
    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( dBdt[i] > 0 )
            dBdt[i] = log10( dBdt[i] );
        else
            dBdt[i] = glob_log_dBdt_min;
    }
    sprintf( cb_label, "(10^x)" );
    sprintf( buf, "./dBdt/dBdt_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, dBdt, 0, para.PicSize, 0, para.PicSize,
            glob_log_dBdt_min, glob_log_dBdt_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "dBdt (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_dBdt_min, glob_log_dBdt_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    fprintf( LogFilefd, "rate of change of magnetic field analysis ... done.\n" );
}



void magnetic_field_analysis() {
    B_analysis();
    divB_analysis();
    dBdt_analysis();
}
