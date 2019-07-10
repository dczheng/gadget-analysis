#include "allvars.h"

void init_img() {
    writelog( "init image ...\n" );
    memset( &image, 0, sizeof( struct image_struct ) );
    mymalloc2( image.props, All.PicSize * sizeof( double ) );
    mymalloc2( image.img, sizeof( double ) * SQR(All.PicSize) );
    mymalloc2( image.num, sizeof( double ) * SQR(All.PicSize) );
    image.img_tmp = image.img;
    image.num_tmp = image.num;
    writelog( "init image ... done.\n" );
    put_sep0;
}

void free_img() {
    myfree( image.props );
    myfree( image.img );
    myfree( image.num );
}

void reset_img() {
    image.img = image.img_tmp;
    image.num = image.num_tmp;
    memset( image.props, 0, All.PicSize * sizeof(double) );
    memset( image.img, 0, SQR(All.PicSize) * sizeof(double) );
    memset( image.num, 0, SQR(All.PicSize) * sizeof(double) );
}

void write_img( char *fn, char *nstr, int mode ) {
    FILE *fd;
    fd = fopen( fn, "w" );
    double v;
    int i, j;
    if ( nstr == NULL ) {
        writelog( "save img ...\n" );
    }
    else {
        writelog( "save `%s` img ...\n", nstr );
    }
    img_min = DBL_MAX;
    img_max = -DBL_MAX;

    for ( i=0; i<All.PicSize2; i++ ) {
        vmin2( img_min, image.img[i], 0 );
        vmax2( img_max, image.img[i]);
    }

    /*
    writelog( "xmin: %g, xmax: %g, ymin: %g, ymax: %g\n",
            img_xmin, img_xmax, img_ymin, img_ymax );

    writelog( "min: %g, max: %g\n", img_min, img_max );
    */

    if ( mode  ) {

        find_global_value( img_min, img_globmin, MPI_DOUBLE, MPI_MIN, MpiComm_Master );
        find_global_value( img_max, img_globmax, MPI_DOUBLE, MPI_MAX, MpiComm_Master );

        find_global_value( img_xmin, img_globxmin, MPI_DOUBLE, MPI_MIN, MpiComm_Master );
        find_global_value( img_xmax, img_globxmax, MPI_DOUBLE, MPI_MAX, MpiComm_Master );
        find_global_value( img_ymin, img_globymin, MPI_DOUBLE, MPI_MIN, MpiComm_Master );
        find_global_value( img_ymax, img_globymax, MPI_DOUBLE, MPI_MAX, MpiComm_Master );

        writelog( "globxmin: %g, globxmax: %g, globymin: %g, glogymax: %g\n",
                    img_globxmin, img_globxmax, img_globymin, img_globymax );

        writelog( "globmin: %g, globmax: %g\n", img_globmin, img_globmax );

    }

    img_proj = All.ProjectDirection;
    img_z = All.RedShift;

    for ( i=0; i<All.PicSize; i++ ) {
        fprintf( fd, "%g ", image.props[i] );
    }
    fprintf( fd, "\n" );

    for ( i=0; i<All.PicSize; i++ ) {
        for ( j=0; j<All.PicSize; j++ ) {
            v = image.img[ i*All.PicSize + j ];
            fprintf( fd, "%g ", v );
            }
        fprintf( fd, "\n" );
    }
    fclose( fd );
}

