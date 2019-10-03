#include "allvars.h"

#ifdef MACHNOISE
void  remove_mach_noise() {
    long i;
    put_header( "remove mach noise" );
    for ( i=0; i<N_Gas; i++ ) {
        if ( i % NTask != ThisTask )
            continue;
        if ( SphP[i].MachNumber > 2 && SphP[i].Density / Time3 / RhoBaryon > 1e5 ) {
            SphP[i].MachNumber = 1;
        }
    }

    do_sync( "" );
    put_end();
}
#endif

