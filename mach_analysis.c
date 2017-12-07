#include "allvars.h"


void mach_analysis() {
    double *mn, mn_max, mn_min, log_mn_max, log_mn_min,
           dx, dy, x, y, glob_log_mn_max, glob_log_mn_min,
           h, dh, lx, ly, bmag, pic_mn_max, pic_mn_min, log_pic_mn_max,
           log_pic_mn_min, glob_log_pic_mn_max, glob_log_pic_mn_min ;
    int i, j, xi, yi, N, Nhalf, i1, i2, j1, j2, li, lj, PicSize;
    char buf[100];
    PicSize = para.PicSize;
    print_log( "mn analysis ..." );
    mn = malloc( sizeof( double ) * PicSize * PicSize );
    memset( mn, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = proj_size / PicSize;
    mn_max = -DBL_MAX;
    mn_min = DBL_MAX;
    pic_mn_max = -DBL_MAX;
    pic_mn_min = DBL_MAX;
    N = para.KernelN;
    Nhalf = N / 2;
    h = para.SofteningTable[0];
    dh = h / Nhalf;
    for ( i=0; i<SliceEnd[0]; i++ ){
        x = P[i].Pos[proj_i];
        y = P[i].Pos[proj_j];
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
                        mn[ i1 * PicSize + j1 ] += SphP[i].MachNumber * KernelMat2D[0][ li*N + lj ] / ( dx * dy );
                }
        }
        else
            mn[ xi * PicSize + yi ] += SphP[i].MachNumber / ( dx * dy );
        if ( SphP[i].MachNumber > mn_max )
            mn_max = SphP[i].MachNumber;
        if ( SphP[i].MachNumber < mn_min )
            mn_min = SphP[i].MachNumber;
    }

    if ( mn_max == -DBL_MAX )
        log_mn_max = -DBL_MAX;
    else
        log_mn_max = log10( mn_max );

    if ( mn_min == DBL_MAX )
        log_mn_min = DBL_MAX;
    else
        log_mn_min = log10( mn_min );

    sprintf( LogBuf, "Particle: mn_max: %g, mn_min: %g\n"
            "log_mn_max: %g, log_mn_min: %g",
            mn_max, mn_min,
            log_mn_max, log_mn_min );
    print_log( LogBuf );

    for ( i=0; i<PicSize*PicSize; i++ ) {
        if ( mn[i] < pic_mn_min && mn[i] > 0 )
            pic_mn_min = mn[i];
        if ( mn[i] > pic_mn_max )
            pic_mn_max = mn[i];
    }

    if ( pic_mn_max == -DBL_MAX )
        log_pic_mn_max = -DBL_MAX;
    else
        log_pic_mn_max = log10( pic_mn_max );

    if ( pic_mn_min == DBL_MAX )
        log_pic_mn_min = DBL_MAX;
    else
        log_pic_mn_min = log10( pic_mn_min );

    sprintf( LogBuf, "Picture: mn_max: %g, mn_min: %g\n"
            "log_mn_max: %g, log_mn_min: %g",
            pic_mn_max, pic_mn_min,
            log_pic_mn_max, log_pic_mn_min );
    print_log( LogBuf );


    if ( ThisTask == 0 )
    if ( access( "./mn/", 0 ) == -1 ){
        print_log( "create directory `./mn` by task 0" );
        if ( mkdir( "./mn", 0755) == -1 ){
            printf( "failed create directory ./mn.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_mn_max, &glob_log_mn_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_mn_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_mn_min, &glob_log_mn_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_mn_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    sprintf( LogBuf, "Particle: glob_log_mn_max: %g, glob_log_mn_min: %g",
            glob_log_mn_max, glob_log_mn_min );
    print_log( LogBuf );
    MPI_Reduce( &log_pic_mn_max, &glob_log_pic_mn_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_pic_mn_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_pic_mn_min, &glob_log_pic_mn_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_pic_mn_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    sprintf( LogBuf, "Picture: glob_log_mn_max: %g, glob_log_mn_min: %g",
            glob_log_pic_mn_max, glob_log_pic_mn_min );
    print_log( LogBuf );

    for ( i=0; i<PicSize*PicSize; i++ ){
        if ( mn[i] > 0 )
            mn[i] = log10( mn[i] );
        else
            mn[i] = glob_log_pic_mn_min;
    }
    sprintf( cb_label, "(10^x) G/Kpc^2" );
    sprintf( buf, "./mn/mn_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( PicSize, PicSize, mn, 0, PicSize, 0, PicSize,
            glob_log_pic_mn_min, glob_log_pic_mn_max, 0, affine );
    sprintf( xlabel, "%g Mpc", proj_size / para.MpcFlag );
    sprintf( title, "mn (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_pic_mn_min, glob_log_pic_mn_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    free( mn );
    print_log( "mn analysis ... done." );
    print_log( sep_str );
}

