#include "allvars.h"

void hg_electrons_analysis() {
    double *hge_n, hge_n_max, hge_n_min, log_hge_n_max, log_hge_n_min,
           dx, dy, x, y, glob_log_hge_n_max, glob_log_hge_n_min,
           h, dh, lx, ly, bmag, pic_hge_n_max, pic_hge_n_min, log_pic_hge_n_max,
           log_pic_hge_n_min, glob_log_pic_hge_n_max, glob_log_pic_hge_n_min, n;
    int i, j, xi, yi, N, Nhalf, i1, i2, j1, j2, li, lj, PicSize;
    char buf[100];
    PicSize = para.PicSize;
    print_log( "hge_n analysis ..." );
    hge_n = malloc( sizeof( double ) * PicSize * PicSize );
    memset( hge_n, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = proj_size / PicSize;
    hge_n_max = -DBL_MAX;
    hge_n_min = DBL_MAX;
    pic_hge_n_max = -DBL_MAX;
    pic_hge_n_min = DBL_MAX;
    N = para.KernelN;
    Nhalf = N / 2;
    h = para.SofteningTable[0];
    dh = h / Nhalf;
    for ( i=0; i<SliceEnd[0]; i++ ){
        x = P[i].Pos[proj_i];
        y = P[i].Pos[proj_j];
        n = SphP[i].CRE_n0 * SphP[i].Density;
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
                        hge_n[ i1 * PicSize + j1 ] += n * KernelMat2D[0][ li*N + lj ] / ( dx * dy );
                }
        }
        else
            hge_n[ xi * PicSize + yi ] += n / ( dx * dy );
        if ( n > hge_n_max && n > 0 )
            hge_n_max = n;
        if ( n < hge_n_min && n > 0 )
            hge_n_min = n;
    }

    if ( hge_n_max == -DBL_MAX )
        log_hge_n_max = -DBL_MAX;
    else
        log_hge_n_max = log10( hge_n_max );

    if ( hge_n_min == DBL_MAX )
        log_hge_n_min = DBL_MAX;
    else
        log_hge_n_min = log10( hge_n_min );

    sprintf( LogBuf, "Particle: hge_n_max: %g, hge_n_min: %g\n"
            "log_hge_n_max: %g, log_hge_n_min: %g",
            hge_n_max, hge_n_min,
            log_hge_n_max, log_hge_n_min );
    print_log( LogBuf );

    for ( i=0; i<PicSize*PicSize; i++ ) {
        if ( hge_n[i] < pic_hge_n_min && hge_n[i] > 0 )
            pic_hge_n_min = hge_n[i];
        if ( hge_n[i] > pic_hge_n_max )
            pic_hge_n_max = hge_n[i];
    }

    if ( pic_hge_n_max == -DBL_MAX )
        log_pic_hge_n_max = -DBL_MAX;
    else
        log_pic_hge_n_max = log10( pic_hge_n_max );

    if ( pic_hge_n_min == DBL_MAX )
        log_pic_hge_n_min = DBL_MAX;
    else
        log_pic_hge_n_min = log10( pic_hge_n_min );

    sprintf( LogBuf, "Picture: hge_n_max: %g, hge_n_min: %g\n"
            "log_hge_n_max: %g, log_hge_n_min: %g",
            pic_hge_n_max, pic_hge_n_min,
            log_pic_hge_n_max, log_pic_hge_n_min );
    print_log( LogBuf );


    if ( ThisTask == 0 )
    if ( access( "./hge_n/", 0 ) == -1 ){
        print_log( "create directory `./hge_n` by task 0" );
        if ( mkdir( "./hge_n", 0755) == -1 ){
            printf( "failed create directory ./hge_n.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_hge_n_max, &glob_log_hge_n_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_hge_n_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_hge_n_min, &glob_log_hge_n_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_hge_n_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    sprintf( LogBuf, "Particle: glob_log_hge_n_max: %g, glob_log_hge_n_min: %g",
            glob_log_hge_n_max, glob_log_hge_n_min );
    print_log( LogBuf );
    MPI_Reduce( &log_pic_hge_n_max, &glob_log_pic_hge_n_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_pic_hge_n_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_pic_hge_n_min, &glob_log_pic_hge_n_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_pic_hge_n_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    sprintf( LogBuf, "Picture: glob_log_hge_n_max: %g, glob_log_hge_n_min: %g",
            glob_log_pic_hge_n_max, glob_log_pic_hge_n_min );
    print_log( LogBuf );

    for ( i=0; i<PicSize*PicSize; i++ ){
        if ( hge_n[i] > 0 )
            hge_n[i] = log10( hge_n[i] );
        else
            hge_n[i] = glob_log_pic_hge_n_min;
    }
    sprintf( cb_label, "(10^x) Kpc^{-2}" );
    sprintf( buf, "./hge_n/hge_n_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( PicSize, PicSize, hge_n, 0, PicSize, 0, PicSize,
            glob_log_pic_hge_n_min, glob_log_pic_hge_n_max, 0, affine );
    sprintf( xlabel, "%g Mpc", proj_size / para.MpcFlag );
    sprintf( title, "hge_n (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_pic_hge_n_min, glob_log_pic_hge_n_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    free( hge_n );
    print_log( "hge_n analysis ... done." );
    print_log( sep_str );
}

