#include "allvars.h"

void write_img( char *fn, int mode ) {
    FILE *fd;
    fd = fopen( fn, "w" );
    double x, y, v, vl1, vl2;
    int i, j, ii, jj, N;
    writelog( "write image data ...\n" );
    N = 0;
    fprintf( fd, "%g ", All.RedShift );
    N += 1;
    fprintf( fd, "%g %g %g %g %g %g ",
            image.xmin, image.xmax,
            image.ymin, image.ymax,
            image.zmin, image.zmax );
    N += 6;

    if ( mode == 0 ){
        fprintf( fd, "%g %g %g %g %g %g %g %g ",
                image.ImgMin, image.ImgMax,
                image.GlobalImgMin, image.GlobalImgMax,
                image.DataMin, image.DataMax,
                image.GlobalDataMin, image.GlobalDataMax );
    }
    else {
        writelog( "log image info: " );
        vl1 = log10( image.GlobalImgMin );
        vl2 = log10( image.GlobalImgMax );
        writelog( "Min=%g, Max=%g\n", vl1, vl2 );
        fprintf( fd, "%g %g %g %g %g %g %g %g ",
                log10(image.ImgMin), log10(image.ImgMax),
                log10(image.GlobalImgMin), log10(image.GlobalImgMax),
                log10(image.DataMin), log10(image.DataMax),
                log10(image.GlobalDataMin), log10(image.GlobalDataMax) );
    }
    N += 8;

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
            if ( mode == 0 )
                fprintf( fd, "%g ", v );
            else {
                if ( v == 0 )
                    fprintf( fd, "%g ", vl1-1 );
                else
                    fprintf( fd, "%g ", log10( v ) );
            }
        }
        fprintf( fd, "\n" );
    }
    fclose( fd );
    writelog( "write image data ... done.\n" );
}

