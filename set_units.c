#include "allvars.h"

void set_units() {
    print_log( "Set Units..." );
    All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
    All.UnitDensity_in_cgs = All.UnitMass_in_g / pow( All.UnitLength_in_cm, 3 );
    All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow( All.UnitLength_in_cm,2 ) / pow( All.UnitTime_in_s, 2 );
    All.G = GRAVITY / pow( All.UnitLength_in_cm, 3 ) * All.UnitMass_in_g * pow( All.UnitTime_in_s, 2 );
    All.Hubble = HUBBLE * All.UnitTime_in_s;

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
    sprintf( LogBuf,  "%-35s: %g", "Gravity constant", All.G );
    sprintf( LogBuf,  "%-35s: %g", "Hubble", All.Hubble );
    print_log( LogBuf );
    if ( All.MpcFlag != 1 ) {
        All.MpcFlag = 1000;
    }
    g2c.cm       = All.UnitLength_in_cm / header.HubbleParam;
    g2c.g        = All.UnitMass_in_g / header.HubbleParam;
    g2c.s        = All.UnitTime_in_s / header.HubbleParam;
    g2c.erg      = g2c.g * SQR(g2c.cm) / SQR( g2c.s );
    print_log( sep_str );
}

