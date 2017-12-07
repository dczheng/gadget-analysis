#include "allvars.h"

void set_units() {
    print_log( "Set Units..." );
    para.UnitTime_in_s = para.UnitLength_in_cm / para.UnitVelocity_in_cm_per_s;
    para.UnitDensity_in_cgs = para.UnitMass_in_g / pow( para.UnitLength_in_cm, 3 );
    para.UnitEnergy_in_cgs = para.UnitMass_in_g * pow( para.UnitLength_in_cm,2 ) / pow( para.UnitTime_in_s, 2 );

    sprintf( LogBuf,  "%-35s: %g", "UnitMass_in_g", para.UnitMass_in_g );
    print_log( LogBuf );
    sprintf( LogBuf,  "%-35s: %g", "UnitTime_in_s", para.UnitTime_in_s );
    print_log( LogBuf );
    sprintf( LogBuf,  "%-35s: %g", "UnitLength_in_cm", para.UnitLength_in_cm );
    print_log( LogBuf );
    sprintf( LogBuf,  "%-35s: %g", "UnitDensity_in_cgs", para.UnitDensity_in_cgs );
    print_log( LogBuf );
    sprintf( LogBuf,  "%-35s: %g", "UnitEnergy_in_cgs", para.UnitEnergy_in_cgs );
    print_log( LogBuf );
    sprintf( LogBuf,  "%-35s: %g", "UnitVelocity_in_cm_per_s", para.UnitVelocity_in_cm_per_s );
    print_log( LogBuf );
    if ( para.MpcFlag != 1 ) {
        para.MpcFlag = 1000;
    }
    print_log( sep_str );
}

