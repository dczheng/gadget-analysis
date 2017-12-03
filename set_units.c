#include "allvars.h"

void set_units() {
    fputs( sep_str, LogFilefd );
    fprintf( LogFilefd, "Set para.Units...\n" );
    para.UnitTime_in_s = para.UnitLength_in_cm / para.UnitVelocity_in_cm_per_s;
    para.UnitDensity_in_cgs = para.UnitMass_in_g / pow( para.UnitLength_in_cm, 3 );
    para.UnitEnergy_in_cgs = para.UnitMass_in_g * pow( para.UnitLength_in_cm,2 ) / pow( para.UnitTime_in_s, 2 );
    fprintf( LogFilefd,  "%-35s: %g\n", "UnitMass_in_g", para.UnitMass_in_g );
    fprintf( LogFilefd,  "%-35s: %g\n", "UnitTime_in_s", para.UnitTime_in_s );
    fprintf( LogFilefd,  "%-35s: %g\n", "UnitLength_in_cm", para.UnitLength_in_cm );
    fprintf( LogFilefd,  "%-35s: %g\n", "UnitDensity_in_cgs", para.UnitDensity_in_cgs );
    fprintf( LogFilefd,  "%-35s: %g\n", "UnitEnergy_in_cgs", para.UnitEnergy_in_cgs );
    fprintf( LogFilefd,  "%-35s: %g\n", "UnitVelocity_in_cm_per_s", para.UnitVelocity_in_cm_per_s );
    if ( para.MpcFlag != 1 ) {
        para.MpcFlag = 1000;
    }
    fputs( sep_str, LogFilefd );
}

