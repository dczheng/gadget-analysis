#include "allvars.h"

void mach_slice() {

    int num, i;
    char buf[100];
    double *data;

    writelog( "gas mach silce ...\n" );
    num = All.SliceEnd[0] - All.SliceStart[0];

    mymalloc2( data, sizeof( double ) * num );

    for ( i=All.SliceStart[0]; i<All.SliceEnd[0]; i++ ) {
        data[i] = SphP[i].MachNumber;
    }

    sprintf( buf, "%sMachNumber", All.OutputDir );
    create_dir( buf );
    sprintf( buf, "%s/MachNumber_%.2f.dat", buf, All.RedShift );

    make_slice_img( 0, data );

    write_img2( buf, "mach number slice" );
    myfree( data );

    writelog( "mach number silce ... done.\n" );
    put_sep;
}
