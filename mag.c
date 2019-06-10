#include "allvars.h"

void mag_slice() {

    int num, i;
    char buf[100];
    double *data;

    writelog( "gas magnetic filed silce ...\n" );
    num = All.SliceEnd[0] - All.SliceStart[0];

    mymalloc2( data, sizeof( double ) * num );

    for ( i=All.SliceStart[0]; i<All.SliceEnd[0]; i++ ) {
        data[i] = get_B( i );
    }

    sprintf( buf, "%sMagneticField", All.OutputDir );
    create_dir( buf );
    sprintf( buf, "%s/MagneticField_%.2f.dat", buf, All.RedShift );

    make_slice_img( 0, data );

    write_img2( buf, "magnetic ield slice" );
    myfree( data );

    writelog( "magnetic field silce ... done.\n" );
    put_sep;
}
