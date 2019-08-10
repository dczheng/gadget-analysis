#include "allvars.h"

void compute_temperature() {
    double yhelium, u, ne, mu, XH;
    int i;

    if ( All.Readu == 0 )
        endruns( "require All.Readu" );
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

