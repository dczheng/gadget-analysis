#include "allvars.h"

void write_img( char *fn ) {
    FILE *fd;
    fd = fopen( fn, "w" );
    double x, y, v, vl1, vl2;
    int i, j, ii, jj, N;
    writelog( "write image data ...\n" );
    N = 0;
    fprintf( fd, "%g ", All.RedShift );
    N += 1;

    fprintf( fd, "%g %g %g %g ",
            image.xmin, image.xmax,
            image.ymin, image.ymax );
    N += 4;

    image.ImgMin = DBL_MAX;
    image.ImgMax = -DBL_MAX;
    for ( i=0; i<All.PicSize2; i++ ) {
        image.ImgMin = vmin( image.img[i], image.ImgMin, 0 );
        image.ImgMax = vmax( image.img[i], image.ImgMax );
    }

    find_global_value( image.ImgMin, image.GlobalImgMin, MPI_DOUBLE, MPI_MIN );
    find_global_value( image.ImgMax, image.GlobalImgMax, MPI_DOUBLE, MPI_MAX );

    fprintf( fd, "%g %g %g %g ",
            image.ImgMin, image.ImgMax,
            image.GlobalImgMin, image.GlobalImgMax );
    N += 4;

    if ( N >= All.PicSize ) {
        printf( "PicSize is too small!\n" );
        endrun( 20180425 );
    }

    for ( i=N; i<All.PicSize; i++ ) {
        fprintf( fd, "0 " );
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

