#include "allvars.h"

void mach_slice() {

    int num, i;
    char buf[100];

    writelog( "gas mach silce ...\n" );
    num = All.SliceEnd[0] - All.SliceStart[0];

    mymalloc2( image.data, sizeof( double ) * num );

    for ( i=All.SliceStart[0]; i<All.SliceEnd[0]; i++ ) {
        image.data[i] = SphP[i].MachNumber;
    }

    sprintf( buf, "%sMachNumber", All.OutputPrefix );
    create_dir( buf );
    sprintf( buf, "%s/MachNumber_%.2f.dat", buf, All.RedShift );

    make_slice_img( 0 );

    write_img2( buf, "mach number slice" );
    myfree( image.data );

    writelog( "mach number silce ... done.\n" );
    put_sep;
}
