#include "allvars.h"

void mag_slice() {

    int num, i;
    char buf[100];
    double dx;

    if ( ThisTask_Local != 0 )
        return;

    writelog( "gas magnetic filed silce ...\n" );
    num = All.SliceEnd[0] - All.SliceStart[0];

    mymalloc2( image.data, sizeof( double ) * num );
    mymalloc2( image.img, sizeof( double ) * SQR( All.PicSize ) );

    for ( i=All.SliceStart[0]; i<All.SliceEnd[0]; i++ ) {
        image.data[i] = get_B( i );
    }

    sprintf( buf, "%sMagneticField", All.OutputPrefix );
    create_dir( buf );
    sprintf( buf, "%s/MagneticField_%.2f.dat", buf, All.RedShift );

    make_slice_img( 0 );

    dx = ( All.End[All.proj_i] - All.Start[All.proj_i] ) / All.PicSize;
    for ( i=0; i<SQR(All.PicSize); i++ ) {
        image.img[i] /= SQR(dx);
    }

    write_img2( buf, "magnetic ield slice" );
    myfree( image.data );
    myfree( image.img );

    writelog( "magnetic field silce ... done.\n" );
    put_sep;
}
