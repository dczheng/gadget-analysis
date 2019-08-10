#include"allvars.h"

void phase() {

    double *phase_x, *phase_y, mm[4];
    long p;
    int flag;

    mymalloc1( phase_x, sizeof(double) * N_Gas );
    mymalloc2( phase_y, sizeof(double) * N_Gas );

    for( p=0; p<N_Gas; p++ ) {
        phase_x[p] = SphP[p].Density / All.RhoBaryon;
        phase_y[p] = SphP[p].Temp;
    }
    writelog( "plot gas phase...\n" );

    flag = 0;
    flag |= PDF2D_BIT_XLOG;
    flag |= PDF2D_BIT_YLOG;

    if ( All.PhaseDensMin > 0 || 
            All.PhaseDensMax > 0 ) {
        if ( All.PhaseDensMin * All.PhaseDensMax == 0 )
            endrun( 20180810 );

        flag |= PDF2D_BIT_FIXEDX;
        mm[0] = All.PhaseDensMin;
        mm[1] = All.PhaseDensMax;

    }

    if ( All.PhaseTempMin > 0 || 
            All.PhaseTempMax > 0 ) {
        if ( All.PhaseTempMin * All.PhaseTempMax == 0 )
            endrun( 20180810 );

        flag |= PDF2D_BIT_FIXEDY;
        mm[2] = All.PhaseTempMin;
        mm[3] = All.PhaseTempMax;

    }

    pdf2d( phase_x, phase_y, NULL, N_Gas, "Phase", flag, mm );
    
    myfree( phase_x );
    myfree( phase_y );
    put_sep;

}
