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
        if ( bmag < B_min && bmag > 0 )
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
    sprintf( xlabel, "%g Mpc", proj_size / para.MpcFlag );
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
    double *divB, divB_max, divB_min, log_divB_max, log_divB_min,
           dx, dy, x, y, glob_log_divB_max, glob_log_divB_min,
           h, dh, lx, ly, bmag, pic_divB_max, pic_divB_min, log_pic_divB_max,
           log_pic_divB_min, glob_log_pic_divB_max, glob_log_pic_divB_min ;
    int i, j, xi, yi, N, Nhalf, i1, i2, j1, j2, li, lj, PicSize;
    char buf[100];
    PicSize = para.PicSize;
    print_log( "div B analysis ..." );
    divB = malloc( sizeof( double ) * PicSize * PicSize );
    memset( divB, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = proj_size / PicSize;
    divB_max = -DBL_MAX;
    divB_min = DBL_MAX;
    pic_divB_max = -DBL_MAX;
    pic_divB_min = DBL_MAX;
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
                        divB[ i1 * PicSize + j1 ] += SphP[i].divB * KernelMat2D[0][ li*N + lj ] / ( dx * dy );
                }
        }
        else
            divB[ xi * PicSize + yi ] += SphP[i].divB / ( dx * dy );
        if ( SphP[i].divB > divB_max )
            divB_max = SphP[i].divB;
        if ( SphP[i].divB < divB_min && SphP[i].divB > 0 )
            divB_min = SphP[i].divB;
    }

    if ( divB_max == -DBL_MAX )
        log_divB_max = -DBL_MAX;
    else
        log_divB_max = log10( divB_max );

    if ( divB_min == DBL_MAX )
        log_divB_min = DBL_MAX;
    else
        log_divB_min = log10( divB_min );

    sprintf( LogBuf, "Particle: divB_max: %g, divB_min: %g\n"
            "log_divB_max: %g, log_divB_min: %g",
            divB_max, divB_min,
            log_divB_max, log_divB_min );
    print_log( LogBuf );

    for ( i=0; i<PicSize*PicSize; i++ ) {
        if ( divB[i] < pic_divB_min && divB[i] > 0 )
            pic_divB_min = divB[i];
        if ( divB[i] > pic_divB_max )
            pic_divB_max = divB[i];
    }

    if ( pic_divB_max == -DBL_MAX )
        log_pic_divB_max = -DBL_MAX;
    else
        log_pic_divB_max = log10( pic_divB_max );

    if ( pic_divB_min == DBL_MAX )
        log_pic_divB_min = DBL_MAX;
    else
        log_pic_divB_min = log10( pic_divB_min );

    sprintf( LogBuf, "Picture: divB_max: %g, divB_min: %g\n"
            "log_divB_max: %g, log_divB_min: %g",
            pic_divB_max, pic_divB_min,
            log_pic_divB_max, log_pic_divB_min );
    print_log( LogBuf );


    if ( ThisTask == 0 )
    if ( access( "./divB/", 0 ) == -1 ){
        print_log( "create directory `./divB` by task 0" );
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
    sprintf( LogBuf, "Particle: glob_log_divB_max: %g, glob_log_divB_min: %g",
            glob_log_divB_max, glob_log_divB_min );
    print_log( LogBuf );
    MPI_Reduce( &log_pic_divB_max, &glob_log_pic_divB_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_pic_divB_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_pic_divB_min, &glob_log_pic_divB_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_pic_divB_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    sprintf( LogBuf, "Picture: glob_log_divB_max: %g, glob_log_divB_min: %g",
            glob_log_pic_divB_max, glob_log_pic_divB_min );
    print_log( LogBuf );

    for ( i=0; i<PicSize*PicSize; i++ ){
        if ( divB[i] > 0 )
            divB[i] = log10( divB[i] );
        else
            divB[i] = glob_log_pic_divB_min;
    }
    sprintf( cb_label, "(10^x) G/Kpc^2" );
    sprintf( buf, "./divB/divB_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( PicSize, PicSize, divB, 0, PicSize, 0, PicSize,
            glob_log_pic_divB_min, glob_log_pic_divB_max, 0, affine );
    sprintf( xlabel, "%g Mpc", proj_size / para.MpcFlag );
    sprintf( title, "divB (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_pic_divB_min, glob_log_pic_divB_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    free( divB );
    print_log( "div B analysis ... done." );
    print_log( sep_str );
}

void dBdt_analysis() {
    double *dBdt, dBdt_max, dBdt_min, log_dBdt_max, log_dBdt_min,
           dx, dy, x, y, glob_log_dBdt_max, glob_log_dBdt_min,
           h, dh, lx, ly, bmag, pic_dBdt_max, pic_dBdt_min, log_pic_dBdt_max,
           log_pic_dBdt_min, glob_log_pic_dBdt_max, glob_log_pic_dBdt_min ;
    int i, j, xi, yi, N, Nhalf, i1, i2, j1, j2, li, lj, PicSize;
    char buf[100];
    PicSize = para.PicSize;
    print_log( "dBdt analysis ..." );
    dBdt = malloc( sizeof( double ) * PicSize * PicSize );
    memset( dBdt, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = proj_size / PicSize;
    dBdt_max = -DBL_MAX;
    dBdt_min = DBL_MAX;
    pic_dBdt_max = -DBL_MAX;
    pic_dBdt_min = DBL_MAX;
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
                        dBdt[ i1 * PicSize + j1 ] += SphP[i].dBdt * KernelMat2D[0][ li*N + lj ] / ( dx * dy );
                }
        }
        else
            dBdt[ xi * PicSize + yi ] += SphP[i].dBdt / ( dx * dy );
        if ( SphP[i].dBdt > dBdt_max )
            dBdt_max = SphP[i].dBdt;
        if ( SphP[i].dBdt < dBdt_min && SphP[i].dBdt > 0 )
            dBdt_min = SphP[i].dBdt;
    }

    if ( dBdt_max == -DBL_MAX )
        log_dBdt_max = -DBL_MAX;
    else
        log_dBdt_max = log10( dBdt_max );

    if ( dBdt_min == DBL_MAX )
        log_dBdt_min = DBL_MAX;
    else
        log_dBdt_min = log10( dBdt_min );

    sprintf( LogBuf, "Particle: dBdt_max: %g, dBdt_min: %g\n"
            "log_dBdt_max: %g, log_dBdt_min: %g",
            dBdt_max, dBdt_min,
            log_dBdt_max, log_dBdt_min );
    print_log( LogBuf );

    for ( i=0; i<PicSize*PicSize; i++ ) {
        if ( dBdt[i] < pic_dBdt_min && dBdt[i] > 0 )
            pic_dBdt_min = dBdt[i];
        if ( dBdt[i] > pic_dBdt_max )
            pic_dBdt_max = dBdt[i];
    }

    if ( pic_dBdt_max == -DBL_MAX )
        log_pic_dBdt_max = -DBL_MAX;
    else
        log_pic_dBdt_max = log10( pic_dBdt_max );

    if ( pic_dBdt_min == DBL_MAX )
        log_pic_dBdt_min = DBL_MAX;
    else
        log_pic_dBdt_min = log10( pic_dBdt_min );

    sprintf( LogBuf, "Picture: dBdt_max: %g, dBdt_min: %g\n"
            "log_dBdt_max: %g, log_dBdt_min: %g",
            pic_dBdt_max, pic_dBdt_min,
            log_pic_dBdt_max, log_pic_dBdt_min );
    print_log( LogBuf );


    if ( ThisTask == 0 )
    if ( access( "./dBdt/", 0 ) == -1 ){
        print_log( "create directory `./dBdt` by task 0" );
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
    sprintf( LogBuf, "Particle: glob_log_dBdt_max: %g, glob_log_dBdt_min: %g",
            glob_log_dBdt_max, glob_log_dBdt_min );
    print_log( LogBuf );
    MPI_Reduce( &log_pic_dBdt_max, &glob_log_pic_dBdt_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_pic_dBdt_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_pic_dBdt_min, &glob_log_pic_dBdt_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_pic_dBdt_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    sprintf( LogBuf, "Picture: glob_log_dBdt_max: %g, glob_log_dBdt_min: %g",
            glob_log_pic_dBdt_max, glob_log_pic_dBdt_min );
    print_log( LogBuf );

    for ( i=0; i<PicSize*PicSize; i++ ){
        if ( dBdt[i] > 0 )
            dBdt[i] = log10( dBdt[i] );
        else
            dBdt[i] = glob_log_pic_dBdt_min;
    }
    sprintf( cb_label, "(10^x) G/Kpc^2" );
    sprintf( buf, "./dBdt/dBdt_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( PicSize, PicSize, dBdt, 0, PicSize, 0, PicSize,
            glob_log_pic_dBdt_min, glob_log_pic_dBdt_max, 0, affine );
    sprintf( xlabel, "%g Mpc", proj_size / para.MpcFlag );
    sprintf( title, "dBdt (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_pic_dBdt_min, glob_log_pic_dBdt_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    free( dBdt );
    print_log( "dBdt analysis ... done." );
    print_log( sep_str );
}

void magnetic_field_analysis() {
    B_analysis();
    divB_analysis();
    dBdt_analysis();
}
