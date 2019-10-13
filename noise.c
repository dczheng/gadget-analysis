#include "allvars.h"

#ifdef MACHNOISE
void  remove_mach_noise() {
    long i;
    put_header( "remove mach noise" );
    for ( i=0; i<N_Gas; i++ ) {
        if ( i % NTask_Local != ThisTask_Local )
            continue;
        if ( SphP[i].MachNumber > 2 && SphP[i].Density / Time3 / RhoBaryon > 1e6 ) {
            SphP[i].MachNumber = 1;
            SphP[i].CRE_C = SphP[i].CRE_Alpha = SphP[i].CRE_qmin = SphP[i].CRE_qmax = 
            SphP[i].CRE_n = SphP[i].CRE_e = 0;
        }
    }

    put_end();
}
#endif

