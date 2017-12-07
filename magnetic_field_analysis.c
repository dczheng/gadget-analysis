#include "allvars.h"

void B_analysis() {
    double *B, B_max, B_min, log_B_max, log_B_min,
           dx, dy, x, y, glob_log_B_max, glob_log_B_min,
           h, dh, lx, ly, bmag, pic_B_max, pic_B_min, log_pic_B_max,
           log_pic_B_min, glob_log_pic_B_max, glob_log_pic_B_min ;
    int i, j, xi, yi, N, Nhalf, i1, i2, j1, j2, li, lj, PicSize;
    char buf[100];
    PicSize = para.PicSize;
    print_log( "magnetic field analysis ..." );
    B = malloc( sizeof( double ) * PicSize * PicSize );
    memset( B, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = proj_size / PicSize;
    B_max = -DBL_MAX;
    B_min = DBL_MAX;
    pic_B_max = -DBL_MAX;
    pic_B_min = DBL_MAX;
    N = para.KernelN;
    Nhalf = N / 2;
    h = para.SofteningTable[0];
    dh = h / Nhalf;
    for ( i=0; i<SliceEnd[0]; i++ ){
        x = P[i].Pos[proj_i];
        y = P[i].Pos[proj_j];
        bmag = sqrt( pow(SphP[i].B[0], 2) + pow(SphP[i].B[1], 2) + pow(SphP[i].B[2], 2.0) );
        //if ( bmag < 1e-5 ) continue;
        if ( bmag > B_max )
            B_max = bmag;
        if ( bmag < B_min )
            B_min = bmag;
        xi = x / dx;
        yi = y / dy;
        i1 = (int)(( x-h ) / dx);
        i2 = (int)(( x+h ) / dx);
        j1 = (int)(( y-h ) / dy);
        j2 = (int)(( y+h ) / dy);
        if ( i1 != xi || i2 != xi || j1 != yi || j2 != yi ) {
            for ( li=0; li<N; li++ )
                for ( lj=0; lj<N; lj++ ){
                        lx = x + ( li-Nhalf ) * dh;
                        ly = y + ( lj-Nhalf ) * dh;
                        i1 = lx / dx;
                        j1 = ly / dy;
                        if ( i1 < 0 || i1 >= PicSize ||
                                j1 < 0 || j1 >= PicSize ) continue;
                        B[ i1 * PicSize + j1 ] += bmag * KernelMat2D[0][ li*N + lj ] / ( dx * dy );
                }
        }
        else
            B[ xi * PicSize + yi ] += bmag / ( dx * dy );
    }

    if ( B_max == -DBL_MAX )
        log_B_max = -DBL_MAX;
    else
        log_B_max = log10( B_max );

    if ( B_min == DBL_MAX )
        log_B_min = DBL_MAX;
    else
        log_B_min = log10( B_min );

    sprintf( LogBuf, "Particle: B_max: %g, B_min: %g\n"
            "log_B_max: %g, log_B_min: %g",
            B_max, B_min,
            log_B_max, log_B_min );
    print_log( LogBuf );

    for ( i=0; i<PicSize*PicSize; i++ ) {
        if ( B[i] < pic_B_min && B[i] > 0 )
            pic_B_min = B[i];
        if ( B[i] > pic_B_max )
            pic_B_max = B[i];
    }

    if ( pic_B_max == -DBL_MAX )
        log_pic_B_max = -DBL_MAX;
    else
        log_pic_B_max = log10( pic_B_max );

    if ( pic_B_min == DBL_MAX )
        log_pic_B_min = DBL_MAX;
    else
        log_pic_B_min = log10( pic_B_min );

    sprintf( LogBuf, "Picture: B_max: %g, B_min: %g\n"
            "log_B_max: %g, log_B_min: %g",
            pic_B_max, pic_B_min,
            log_pic_B_max, log_pic_B_min );
    print_log( LogBuf );


    if ( ThisTask == 0 )
    if ( access( "./B/", 0 ) == -1 ){
        print_log( "create directory `./B` by task 0" );
        if ( mkdir( "./B", 0755) == -1 ){
            printf( "failed create directory ./B.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_B_max, &glob_log_B_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_B_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_B_min, &glob_log_B_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_B_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    sprintf( LogBuf, "Particle: glob_log_B_max: %g, glob_log_B_min: %g",
            glob_log_B_max, glob_log_B_min );
    print_log( LogBuf );
    MPI_Reduce( &log_pic_B_max, &glob_log_pic_B_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_pic_B_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_pic_B_min, &glob_log_pic_B_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_pic_B_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    sprintf( LogBuf, "Picture: glob_log_B_max: %g, glob_log_B_min: %g",
            glob_log_pic_B_max, glob_log_pic_B_min );
    print_log( LogBuf );

    for ( i=0; i<PicSize*PicSize; i++ ){
        if ( B[i] > 0 )
            B[i] = log10( B[i] );
        else
            B[i] = glob_log_pic_B_min;
    }
    sprintf( cb_label, "(10^x) G/Kpc^2" );
    sprintf( buf, "./B/B_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( PicSize, PicSize, B, 0, PicSize, 0, PicSize,
            glob_log_pic_B_min, glob_log_pic_B_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "B (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_pic_B_min, glob_log_pic_B_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    free( B );
    print_log( "magnetic field analysis ... done." );
    print_log( sep_str );
}

void divB_analysis() {

}

void dBdt_analysis() {

}


void magnetic_field_analysis() {
    B_analysis();
    divB_analysis();
    dBdt_analysis();
}
