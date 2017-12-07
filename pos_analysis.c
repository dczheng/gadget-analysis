#include "allvars.h"

void pos_analysis( int pt ){
    FILE *fd;
    int i,j, xi, yi;
    double *rho, x, y, dx, dy, rho_max, rho_min,
           log_rho_max, log_rho_min,
           glob_log_rho_max, glob_log_rho_min;
    char buf[200], buf1[100];
    long num, offset;
    sprintf( LogBuf, "particle %d positin analysis ...", pt );
    print_log( LogBuf );
    rho = ( double* ) malloc( sizeof(double) * para.PicSize * para.PicSize );
    memset( rho, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
    num = header.npartTotal[pt];
    rho_max = -DBL_MAX;
    rho_min = DBL_MAX;
    num += ( (long long)header.npartTotalHighWord[pt] ) << 32;
    offset = 0;
    for ( i=0; i<pt; i++ ) {
        offset = header.npartTotal[i];
        offset += ( (long long)header.npartTotalHighWord[i] ) << 32;
    }
    for ( i=0; i<num; i++ ) {
        x = P[offset+i].Pos[0];
        y = P[offset+i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        if ( header.mass[pt] == 0 )
            rho[ xi*para.PicSize + yi ] += P[offset+i].Mass / (dx*dy*BoxSize) / ( cgs_g/CUBE(cgs_cm) );
        else
            rho[ xi*para.PicSize + yi ] += header.mass[pt] / (dx*dy*BoxSize) / ( cgs_g/CUBE(cgs_cm) );
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

    sprintf( LogBuf, "rho_max: %g, rho_min: %g\n"
            "log_rho_max: %g, log_rho_min: %g",
            rho_max, rho_min, log_rho_max, log_rho_min );
    print_log( LogBuf );
    sprintf( buf1, "./rho_%d", pt );
    if ( ThisTask == 0 )
    if ( access( buf1, 0 ) == -1 ){
        sprintf( LogBuf, "create directory `%s` by task %d", buf1, ThisTask );
        print_log( LogBuf );
        if ( mkdir( buf1, 0755) == -1 ){
            printf( "failed create directory %s.\n", buf1 );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_rho_max, &glob_log_rho_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_rho_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_rho_min, &glob_log_rho_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_rho_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    sprintf( LogBuf, "glob_log_rho_max: %g, glob_log_rho_min: %g",
            glob_log_rho_max, glob_log_rho_min );
    print_log( LogBuf );
    for ( i=0; i<para.PicSize*para.PicSize; i++ )
        if ( rho[i] > 0 )
            rho[i] = log10( rho[i] );
        else
            rho[i] = glob_log_rho_min;

    sprintf( cb_label, "(10^x)   g/cm^3" );
    sprintf( buf, "%s/rho_%d_%.2f.png", buf1, pt, RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, rho, 0, para.PicSize, 0, para.PicSize, glob_log_rho_min, glob_log_rho_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "particle %d density (z=%.2f)", pt, RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_rho_min, glob_log_rho_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;


    free( rho );
    sprintf( LogBuf, "particle %d positin analysis ... done.", pt );
    print_log( LogBuf );
    print_log( sep_str );
}

