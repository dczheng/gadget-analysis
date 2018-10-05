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
    put_block_line;

}

void gas_temperature() {

    int num, i, PicSize2;
    char buf[100];
    writelog( "gas temperature silce ...\n" );
    PicSize2 = All.PicSize2;
    num = All.SliceEnd[0] - All.SliceStart[0];
    mymalloc2( image.data, sizeof( double ) * num );
    mymalloc2( image.img, sizeof( double ) * SQR( All.PicSize ) );

    for ( i=All.SliceStart[0]; i<num; i++ ) {
        image.data[i] = SphP[i].Temp;
    }

    create_dir( "./gas_temperature" );
    sprintf( buf, "./gas_temperature/%s_gas_temperature_%03i_%.2f.dat", All.FilePrefix, ThisTask, All.RedShift );

    make_slice_img( 0 );

    for ( i=0; i<PicSize2; i++ ){
        image.img[i] *= All.UnitLength_in_cm;
    }

    write_img2( buf, "gas temperature slice" );

    myfree( image.data );
    myfree( image.img );
    writelog( "gas Temperature silce ... done.\n" );
    put_block_line;

}
