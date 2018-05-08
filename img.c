#include "allvars.h"

void init_img() {
    writelog( "init image ...\n" );
    memset( &image, 0, sizeof( struct image_struct ) );
    mymalloc( image.props, All.PicSize * sizeof( double ) );
    writelog( "init image ... done.\n" );
}

void free_img() {
    myfree( image.props );
}

void write_img( char *fn ) {
    FILE *fd;
    fd = fopen( fn, "w" );
    double v;
    int i, j;
    writelog( "write image data ...\n" );
    img_min = DBL_MAX;
    img_max = -DBL_MAX;
    for ( i=0; i<All.PicSize2; i++ ) {
        img_min = vmin( image.img[i], img_min, 0 );
        img_max = vmax( image.img[i], img_max );
    }

    find_global_value( img_min, img_globmin, MPI_DOUBLE, MPI_MIN );
    find_global_value( img_max, img_globmax, MPI_DOUBLE, MPI_MAX );

    writelog( "min: %g, max: %g, globmin: %g, globmax: %g\n",
            img_min, img_max, img_globmin, img_globmax );

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
    writelog( "write image data ... done.\n" );
}

