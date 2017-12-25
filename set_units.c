#include "allvars.h"

void set_units() {
    print_log( "Set Units..." );
    All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
    All.UnitDensity_in_cgs = All.UnitMass_in_g / pow( All.UnitLength_in_cm, 3 );
    All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow( All.UnitLength_in_cm,2 ) / pow( All.UnitTime_in_s, 2 );

    sprintf( LogBuf,  "%-35s: %g", "UnitMass_in_g", All.UnitMass_in_g );
    print_log( LogBuf );
    sprintf( LogBuf,  "%-35s: %g", "UnitTime_in_s", All.UnitTime_in_s );
    print_log( LogBuf );
    sprintf( LogBuf,  "%-35s: %g", "UnitLength_in_cm", All.UnitLength_in_cm );
    print_log( LogBuf );
    sprintf( LogBuf,  "%-35s: %g", "UnitDensity_in_cgs", All.UnitDensity_in_cgs );
    print_log( LogBuf );
    sprintf( LogBuf,  "%-35s: %g", "UnitEnergy_in_cgs", All.UnitEnergy_in_cgs );
    print_log( LogBuf );
    sprintf( LogBuf,  "%-35s: %g", "UnitVelocity_in_cm_per_s", All.UnitVelocity_in_cm_per_s );
    print_log( LogBuf );
    if ( All.MpcFlag != 1 ) {
        All.MpcFlag = 1000;
    }
    print_log( sep_str );
}

