#include"allvars.h"

void phase() {
#ifdef PHASE

    double *phase_x, *phase_y, mm[4];
    long p, index;
    int flag;

    put_header( "gas phase" );

    mymalloc1( phase_x, sizeof(double) * N_Gas );
    mymalloc2( phase_y, sizeof(double) * N_Gas );

    for( p=0, index=0; p<N_Gas; p++ ) {
#ifdef PHASENOSFR
        if ( SphP[p].sfr>0 )
            continue;
#endif
        phase_x[index] = SphP[p].Density / CUBE(Time) / RhoBaryon;
        phase_y[index] = SphP[p].Temp;
        index ++;
    }

    flag = 0;
    flag |= PDF2D_BIT_XLOG;
    flag |= PDF2D_BIT_YLOG;
    //flag |= PDF2D_BIT_UNIT_AREA;

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

    pdf2d( phase_x, phase_y, NULL, index, "Phase", flag, mm );
    
    myfree( phase_x );
    myfree( phase_y );

#endif
}
