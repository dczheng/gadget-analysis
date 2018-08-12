#include "allvars.h"

void init_img() {
    writelog( "init image ...\n" );
    memset( &image, 0, sizeof( struct image_struct ) );
    mymalloc2( image.props, All.PicSize * sizeof( double ) );
    writelog( "init image ... done.\n" );
    put_block_line;
}

void free_img() {
    myfree( image.props );
}

void write_img( char *fn, char *s, int mode ) {
    FILE *fd;
    fd = fopen( fn, "w" );
    double v;
    int i, j;
    writelog( "write `%s` data ...\n", s );
    img_min = DBL_MAX;
    img_max = -DBL_MAX;

    for ( i=0; i<All.PicSize2; i++ ) {
        img_min = vmin( image.img[i], img_min, 0 );
        img_max = vmax( image.img[i], img_max );
    }

    writelog( "xmin: %g, xmax: %g, ymin: %g, ymax: %g\n",
            img_xmin, img_xmax, img_ymin, img_ymax );

    writelog( "min: %g, max: %g\n", img_min, img_max );

    if ( mode  ) {

        find_global_value( img_min, img_globmin, MPI_DOUBLE, MPI_MIN );
        find_global_value( img_max, img_globmax, MPI_DOUBLE, MPI_MAX );

        find_global_value( img_xmin, img_globxmin, MPI_DOUBLE, MPI_MIN );
        find_global_value( img_xmax, img_globxmax, MPI_DOUBLE, MPI_MAX );
        find_global_value( img_ymin, img_globymin, MPI_DOUBLE, MPI_MIN );
        find_global_value( img_ymax, img_globymax, MPI_DOUBLE, MPI_MAX );

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

