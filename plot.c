#include "allvars.h"

int cpn;
double *cp, *red, *green, *blue, affine[6];
char cb_s;
FILE *fp_tmp;

void show_plot_info() {
    sprintf( LogBuf, "plot info:\n%-25s: %i\n"
                                 "%-25s: %i\n"
                                 "%-25s: %i\n"
                                 "%-25s: %i\n"
                                 "%-25s: %s\n"
                                 "%-25s: %s\n"
                                 "%-25s: %s\n"
                                 "%-25s: %s\n"
                                 "%-25s: %s\n"
                                 "%-25s: %g",
                                 "log_flag", plot_info.log_flag,
                                 "global_colorbar_flag", plot_info.global_colorbar_flag,
                                 "istart", plot_info.istart,
                                 "iend", plot_info.iend,
                                 "data_name", plot_info.data_name,
                                 "cb_label", plot_info.cb_label,
                                 "xlabel", plot_info.xlabel,
                                 "ylabel", plot_info.ylabel,
                                 "title", plot_info.title,
                                 "h", plot_info.h );
    print_log( LogBuf );
    print_log( sep_str );
}

void init_plot() {
    int i;
    print_log( "initialize plot..." );
    affine[0] = 1;
    affine[1] = 0;
    affine[2] = 0;
    affine[3] = 1;
    affine[4] = 0;
    affine[5] = 0;
    cb_s = 'R';
    cpn = 30;
    cp = malloc( sizeof(double) * cpn );
    red = malloc( sizeof(double) * cpn );
    green = malloc( sizeof(double) * cpn );
    blue = malloc( sizeof(double) * cpn );
    for ( i=0; i<cpn; i++ ) {
        cp[i] = i / (double)cpn;
        red[i] =   pow( i, 0.9 ) / pow( cpn, 0.9 );
        green[i] =  pow( i, 2 ) / pow( cpn, 2 );
        blue[i] =    pow( i, 0.5 ) / pow( cpn, 0.5 );
    }

    plot_info.log_flag = 1;
    plot_info.global_colorbar_flag = 1;
    plot_info.istart = 0;
    plot_info.iend = 0;
    plot_info.h = 0;
    sprintf( plot_info.data_name, "none" );
    sprintf( plot_info.cb_label, "" );
    sprintf( plot_info.xlabel, "" );
    sprintf( plot_info.ylabel, "" );
    sprintf( plot_info.title, "" );
    show_plot_info();
    print_log( sep_str );
}

void plot_slice() {
    double *img,
           data_max, data_min, log_data_max, log_data_min,
           glob_data_max, glob_data_min, glob_log_data_max, glob_log_data_min,
           img_max, img_min, log_img_max, log_img_min,
           glob_img_max, glob_img_min, glob_log_img_max, glob_log_img_min,
           dx, dy, x, y, h, dh, lx, ly, v,
           cb_max, cb_log_max, cb_min, cb_log_min;
    int i, j, xi, yi, N, Nhalf, i1, i2, j1, j2, li, lj, PicSize;
    char buf[100];
    PicSize = All.PicSize;
    sprintf( LogBuf, "plot %s slice ...", plot_info.data_name );
    print_log( LogBuf );
    img = malloc( sizeof( double ) * PicSize * PicSize );
    memset( img, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = proj_size / PicSize;
    data_max = DBL_MIN;
    data_min = DBL_MAX;
    img_max = DBL_MIN;
    img_min = DBL_MAX;
    N = All.KernelN;
    Nhalf = N / 2;
    h = plot_info.h;
    dh = h / Nhalf;
    for ( i=plot_info.istart; i<plot_info.iend; i++ ){
        x = P[i].Pos[proj_i];
        y = P[i].Pos[proj_j];
        x -= All.Start[proj_i];
        y -= All.Start[proj_j];
        //printf( "%g %g %g\n", x, y, proj_size );
        v = plot_info.data[ i-plot_info.istart ];
        xi = x / dx;
        yi = y / dy;
        xi = ( xi >= PicSize-1 ) ? PicSize-1 : xi;
        yi = ( yi >= PicSize-1 ) ? PicSize-1 : yi;
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
                        img[ i1 * PicSize + j1 ] += v * KernelMat2D[0][ li*N + lj ] / ( dx * dy );
                }
        }
        else
            img[ xi * PicSize + yi ] += v / ( dx * dy );

        if ( plot_info.data[i-plot_info.istart] > 0 ){
            data_max = ( v > data_max ) ? v : data_max;
            data_min = ( v < data_min ) ? v : data_min;
        }

    }

    for ( i=0; i<PicSize*PicSize; i++ )
        if ( img[i] > 0 ){
            img_max = ( img[i] > img_max ) ? img[i] : img_max;
            img_min = ( img[i] < img_min ) ? img[i] : img_min;
        }

    log_data_max = log10( data_max );
    log_data_min = log10( data_min );
    log_img_max = log10( img_max );
    log_img_min = log10( img_min );

    MPI_Reduce( &data_max, &glob_data_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &data_min, &glob_data_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_data_max, &glob_log_data_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_data_min, &glob_log_data_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );

    MPI_Reduce( &img_max, &glob_img_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &img_min, &glob_img_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_img_max, &glob_log_img_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_img_min, &glob_log_img_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );

    MPI_Barrier( MPI_COMM_WORLD );

    MPI_Bcast( &glob_img_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_img_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_img_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_img_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    MPI_Bcast( &glob_img_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_img_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_img_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_img_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    MPI_Barrier( MPI_COMM_WORLD );

    print_log( sep_str );
    sprintf( LogBuf, "Particle info: \nmax: %g, min: %g, log_max: %g, log_min: %g\n"
            "Image info: \nmax: %g, min: %g, log_max: %g, log_min: %g",
            data_max, data_min, log_data_max, log_data_min,
            img_max, img_min, log_img_max, log_img_min );
    print_log( LogBuf );
    print_log( sep_str );
    sprintf( LogBuf, "Global particle info: \nmax: %g, min: %g, log_max: %g, log_min: %g\n"
            "Global image info: \nmax: %g, min: %g, log_max: %g, log_min: %g",
            glob_data_max, glob_data_min, glob_log_data_max, glob_log_data_min,
            glob_img_max, glob_img_min, glob_log_img_max, glob_log_img_min );
    print_log( LogBuf );
    print_log( sep_str );

    if ( ThisTask == 0 ){
        sprintf( buf, "./%s/", plot_info.data_name );
        if ( access( buf, 0 ) == -1 ){
            sprintf( LogBuf, "create directory `./%s/` by task 0", plot_info.data_name );
            print_log( LogBuf );
            if ( mkdir( buf, 0755) == -1 ){
                printf( "failed create directory %s.\n", buf );
                endrun( 20171130 );
            }
        }
    }

    if ( plot_info.log_flag == 1 ){
        for ( i=0; i<PicSize*PicSize; i++ ){
            if ( img[i] > 0 )
                img[i] = log10( img[i] );
            else
             img[i] = glob_log_img_min;
        }
    }

    sprintf( buf, "./%s/%s_%.2f\n", plot_info.data_name, plot_info.data_name, RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    if ( plot_info.global_colorbar_flag == 1 ){
        cb_max = glob_img_max;
        cb_min = glob_img_min;
        cb_log_max = glob_log_img_max;
        cb_log_min = glob_log_img_min;
    }
    else{
        cb_max = img_max;
        cb_min = img_min;
        cb_log_max = log_img_max;
        cb_log_min = log_img_min;
    }
    if ( plot_info.log_flag == 1 ){
        giza_render( PicSize, PicSize, img, 0, PicSize, 0, PicSize,
            glob_log_img_min, glob_log_img_max, 0, affine );
        giza_colour_bar( &cb_s, 1, 3, cb_log_min, cb_log_max, plot_info.cb_label );
    }
    else{
        giza_render( PicSize, PicSize, img, 0, PicSize, 0, PicSize,
            glob_img_min, glob_img_max, 0, affine );
        giza_colour_bar( &cb_s, 1, 3, cb_min, cb_max, plot_info.cb_label );
    }
    if ( strcmp( plot_info.title, "" ) == 0 ){
        sprintf( plot_info.title, "%s (z=%.2f)", plot_info.data_name, RedShift );
    }
    giza_label( plot_info.xlabel, plot_info.ylabel, plot_info.title );
    show_plot_info();

    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    free( img );
    sprintf( LogBuf, "plot %s slice ... done", plot_info.data_name );
    print_log( LogBuf );
    print_log( sep_str );
}


void free_plot() {
    print_log( "free plot ..." );
    free( cp );
    free( red );
    free( green );
    free( blue );
}
