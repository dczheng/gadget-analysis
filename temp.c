#include "allvars.h"

void compute_temperature() {
    double yhelium, u, ne, mu, XH;
    int i;

    writelog( "compute gas temprature...\n" );
    XH = HYDROGEN_MASSFRAC;
    yhelium = ( 1 - XH ) / ( 4 * XH );
    for ( i=0; i<N_Gas; i++ ) {
        u = SphP[i].u * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
        ne = SphP[i].elec;
        mu = ( 1 + 4 * yhelium ) / ( 1 + yhelium + ne );
        SphP[i].Temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
    }
    writelog( "compute gas temprature... done.\n" );
    put_sep;

}

void temperature_slice() {

    int num, i, PicSize2;
    char buf[100];
    double *data;

    writelog( "gas temperature silce ...\n" );
    PicSize2 = All.PicSize2;
    num = All.SliceEnd[0] - All.SliceStart[0];

    mymalloc2( data, sizeof( double ) * num );

    for ( i=All.SliceStart[0]; i<num; i++ ) {
        data[i] = SphP[i].Temp;
    }

    sprintf( buf, "%sTemp", All.OutputDir );
    create_dir( buf );
    sprintf( buf, "%s/Temp_%.2f.dat", buf, All.RedShift );

    make_slice_img( 0, data );

    for ( i=0; i<PicSize2; i++ ){
        image.img[i] *= All.UnitLength_in_cm;
    }

    write_img2( buf, "gas temperature slice" );

    myfree( data );
    writelog( "gas Temperature silce ... done.\n" );
    put_sep;

}
