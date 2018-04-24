#include "allvars.h"

int cpn;
double *cp, *red, *green, *blue, affine[6];
char cb_s;
FILE *fp_tmp;

void show_pli() {
    writelog( "plot info:\n%-25s: %i\n"
                                 "%-25s: %i\n"
                                 "%-25s: %i\n"
                                 "%-25s: %i\n"
                                 "%-25s: %i\n"
                                 "%-25s: %i\n"
                                 "%-25s: %s\n"
                                 "%-25s: %s\n"
                                 "%-25s: %s\n"
                                 "%-25s: %s\n"
                                 "%-25s: %s\n"
                                 "%-25s: %g\n",
                                 "log_flag", pli.log_flag,
                                 "tick_flag", pli.tick_flag,
                                 "global_colorbar_flag", pli.global_colorbar_flag,
                                 "istart", pli.istart,
                                 "iend", pli.iend,
                                 "PicSize", pli.PicSize,
                                 "data_name", pli.data_name,
                                 "cb_label", pli.cb_label,
                                 "xlabel", pli.xlabel,
                                 "ylabel", pli.ylabel,
                                 "title", pli.title,
                                 "h", pli.h );
    writelog( sep_str );
}

void init_plot() {
    int i;
    writelog( "initialize plot...\n" );
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

    pli.log_flag = 0;
    pli.tick_flag = 0;
    pli.global_colorbar_flag = 1;
    pli.istart = 0;
    pli.iend = 0;
    pli.h = 0;
    pli.PicSize = All.PicSize;
    sprintf( pli.data_name, "none" );
    sprintf( pli.cb_label, "" );
    sprintf( pli.xlabel, "" );
    sprintf( pli.ylabel, "" );
    sprintf( pli.title, "" );

    show_pli();
    writelog( sep_str );
}

void plot_imshow() {

    int i, PicSize;
    double *img, ImgMax, ImgMin, LogImgMax, LogImgMin,
           GlobImgMax, GlobImgMin, GlobLogImgMin, GlobLogImgMax,
           cb_max, cb_min, cb_log_min, cb_log_max;
    FILE *fp_tmp;
    char buf[500];

    ImgMax = DBL_MIN;
    ImgMin = DBL_MAX;

    PicSize = pli.PicSize;
    img = pli.img;

    for ( i=0; i<SQR(PicSize); i++ ) {
        ImgMax = ( img[i] > ImgMax ) ? img[i] : ImgMax;
        ImgMin = ( (img[i] < ImgMin && img[i] > 0 )) ? img[i] : ImgMin;
    }
    LogImgMax = log10( ImgMax );
    LogImgMin = log10( ImgMin );

    MPI_Reduce( &ImgMax, &GlobImgMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &ImgMin, &GlobImgMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &LogImgMax, &GlobLogImgMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &LogImgMin, &GlobLogImgMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );

    writelog( sep_str );
    writelog( "Image info: \nLocal: max: %g, min: %g, log_max: %g, log_min: %g\n",
            ImgMax, ImgMin, LogImgMax, LogImgMin );
    writelog( "Global: max: %g, min: %g, log_max: %g, log_min: %g\n",
            GlobImgMax, GlobImgMin, GlobLogImgMax, GlobLogImgMin );
    writelog( sep_str );

    if ( ThisTask == 0 ){
        sprintf( buf, "./%s/", pli.data_name );
        if ( access( buf, 0 ) == -1 ){
            writelog( "create directory `./%s/` by task 0\n", pli.data_name );
            if ( mkdir( buf, 0755) == -1 ){
                printf( "failed create directory %s.\n", buf );
                endrun( 20180102 );
            }
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );

    sprintf( buf, "./%s/%s_%.2f", pli.data_name, pli.data_name, All.RedShift );

    if ( pli.global_colorbar_flag == 1 ){
        cb_max = GlobImgMax;
        cb_min = GlobImgMin;
        cb_log_max = GlobLogImgMax;
        cb_log_min = GlobLogImgMin;
    }
    else{
        cb_max = ImgMax;
        cb_min = ImgMin;
        cb_log_max = LogImgMax;
        cb_log_min = LogImgMin;
    }

    if ( pli.log_flag == 1 )
        for ( i=0; i<SQR(PicSize); i++ ) {
            if ( img[i] > 0 )
                img[i] = log10( img[i] );
            else
                img[i] = LogImgMin;
        }

    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );

    if ( pli.log_flag == 1 ){
        giza_render( PicSize, PicSize, img, 0, PicSize, 0, PicSize,
            GlobLogImgMin, GlobLogImgMax, 0, affine );
        giza_colour_bar( &cb_s, 1, 3, cb_log_min, cb_log_max, pli.cb_label );
    }
    else{
        giza_render( PicSize, PicSize, img, 0, PicSize, 0, PicSize,
            GlobImgMin, GlobImgMax, 0, affine );
        giza_colour_bar( &cb_s, 1, 3, cb_min, cb_max, pli.cb_label );
    }

    if ( strcmp( pli.title, "" ) == 0 ){
        sprintf( pli.title, "%s (z=%.2f)", pli.data_name, All.RedShift );
    }

    giza_label( pli.xlabel, pli.ylabel, pli.title );

    if ( pli.tick_flag ){
        giza_box( &pli.box_xopt,
                   pli.box_xtick,
                   pli.box_ytick,
                  &pli.box_yopt,
                   pli.box_ytick,
                   pli.box_nysub );
    }

    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;

    show_pli();

}

void plot_slice() {
    double *img,
           DataMax, DataMin, LogDataMax, LogDataMin,
           GlobDataMax, GlobDataMin, GlobLogDataMax, GlobLogDataMin,
           dx, dy, x, y, h, dh, lx, ly, v;
    int i, j, xi, yi, N, Nhalf, i1, i2, j1, j2, li, lj, PicSize;
    PicSize = pli.PicSize;
    writelog( "plot %s slice ...\n", pli.data_name );
    img = pli.img;
    memset( img, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = proj_size / PicSize;
    DataMax = DBL_MIN;
    DataMin = DBL_MAX;
    N = All.KernelN;
    Nhalf = N / 2;
    h = pli.h;
    dh = h / Nhalf;
    for ( i=pli.istart; i<pli.iend; i++ ){
        x = P[i].Pos[proj_i];
        y = P[i].Pos[proj_j];
        x -= All.Start[proj_i];
        y -= All.Start[proj_j];
        //printf( "%g %g %g\n", x, y, proj_size );
        v = pli.data[ i-pli.istart ];
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

        if ( pli.data[i-pli.istart] > 0 ){
            DataMax = ( v > DataMax ) ? v : DataMax;
            DataMin = ( v < DataMin && v > 0 ) ? v : DataMin;
        }

    }

    LogDataMax = log10( DataMax );
    LogDataMin = log10( DataMin );

    MPI_Reduce( &DataMax, &GlobDataMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &DataMin, &GlobDataMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &LogDataMax, &GlobLogDataMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &LogDataMin, &GlobLogDataMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );

    writelog( sep_str );
    writelog( "Particle info: \nLocal: max: %g, min: %g, log_max: %g, log_min: %g\n",
            DataMax, DataMin, LogDataMax, LogDataMin );
    writelog( "Global: max: %g, min: %g, log_max: %g, log_min: %g\n",
            GlobDataMax, GlobDataMin, GlobLogDataMax, GlobLogDataMin );
    writelog( sep_str );

    plot_imshow();

    writelog( "plot %s slice ... done\n", pli.data_name );
    writelog( sep_str );

}

void free_plot() {
    writelog( "free plot ...\n" );
    free( cp );
    free( red );
    free( green );
    free( blue );
}
