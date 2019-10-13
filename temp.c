#include "allvars.h"

#ifdef COMPUTETEMP
void compute_temperature() {

#ifdef READTEMP
        return;
#else

    double yhelium, u, ne, mu, XH;
    long i;

    writelog( "compute gas temprature...\n" );
    XH = HYDROGEN_MASSFRAC;
    yhelium = ( 1 - XH ) / ( 4 * XH );
    for ( i=0; i<N_Gas; i++ ) {

        if ( i % NTask_Local != ThisTask_Local )
            continue;

        u = SphP[i].u * UnitEnergy_in_cgs / All.UnitMass_in_g;
        //u = u * 0.5;
#ifndef READELEC
            ne = 1 + 2 * yhelium;
#else
            ne = SphP[i].elec;
#endif

        mu = ( 1 + 4 * yhelium ) / ( 1 + yhelium + ne );
        SphP[i].Temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

    }
    writelog( "compute gas temprature... done.\n" );
#endif

}

#endif
