#include "allvars.h"

void compute_temperature() {
    double yhelium, u, ne, mu, XH;
    long i;

    if ( All.ReadTemp ) {
        return;
    }
    if ( All.Readu == 0  )
        endruns( "require All.Readu" );
    writelog( "compute gas temprature...\n" );
    XH = HYDROGEN_MASSFRAC;
    yhelium = ( 1 - XH ) / ( 4 * XH );
    for ( i=0; i<N_Gas; i++ ) {

        u = SphP[i].u * UnitEnergy_in_cgs / All.UnitMass_in_g;
        //u = u * 0.5;
        if ( All.ReadElec == 0  )
            ne = 1 + 2 * yhelium;
        else
            ne = SphP[i].elec;

        mu = ( 1 + 4 * yhelium ) / ( 1 + yhelium + ne );
        SphP[i].Temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

    }
    writelog( "compute gas temprature... done.\n" );
    put_sep;

}

