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
}

void free_img() {
    reset_img();
    myfree( image.props );
    myfree( image.img );
    myfree( image.num );
}

void reset_img() {

    int i;
    image.img = image.img_tmp;
    image.num = image.num_tmp;
    memset( image.props, 0, All.PicSize * sizeof(double) );
    memset( image.img, 0, SQR(All.PicSize) * sizeof(double) );
    memset( image.num, 0, SQR(All.PicSize) * sizeof(double) );

    for( i=0; i<IMG_PROPS_START; i++ )
        image.props[i] = 0;

}

//#define WRITE_IMG_DEBUG
void write_img( char *fn ) {

    FILE *fd;
    fd = myfopen( "w", fn );
    double v;
    int i, j;

    img_proj = All.ProjectDirection;
    img_z = Redshift;

    for ( i=0; i<All.PicSize; i++ ) {
        fprintf( fd, "%g ", image.props[i] );
#ifdef WRITE_IMG_DEBUG
        printf( "%g ", image.props[i] );
#endif
    }
    fprintf( fd, "\n" );
#ifdef WRITE_IMG_DEBUG
    printf( "\n" );
    endruns( "write_img-debuf" );
#endif

    for ( i=0; i<All.PicSize; i++ ) {
        for ( j=0; j<All.PicSize; j++ ) {
            v = image.img[ i*All.PicSize + j ];
            fprintf( fd, "%e ", v );
            }
        fprintf( fd, "\n" );
    }
    fclose( fd );
}

