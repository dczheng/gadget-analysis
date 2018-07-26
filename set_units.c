#include "allvars.h"

void set_units() {
    writelog( "Set Units...\n" );
    All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
    All.UnitDensity_in_cgs = All.UnitMass_in_g / pow( All.UnitLength_in_cm, 3 );
    All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow( All.UnitLength_in_cm,2 ) / pow( All.UnitTime_in_s, 2 );
    All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow( All.UnitTime_in_s, 2 );
    All.G = GRAVITY / pow( All.UnitLength_in_cm, 3 ) * All.UnitMass_in_g * pow( All.UnitTime_in_s, 2 );
    All.Hubble = HUBBLE * All.UnitTime_in_s;

    writelog( "%-35s: %g\n", "UnitMass_in_g", All.UnitMass_in_g );
    writelog( "%-35s: %g\n", "UnitTime_in_s", All.UnitTime_in_s );
    writelog( "%-35s: %g\n", "UnitLength_in_cm", All.UnitLength_in_cm );
    writelog( "%-35s: %g\n", "UnitDensity_in_cgs", All.UnitDensity_in_cgs );
    writelog( "%-35s: %g\n", "UnitEnergy_in_cgs", All.UnitEnergy_in_cgs );
    writelog( "%-35s: %g\n", "UnitPressure_in_cgs", All.UnitPressure_in_cgs );
    writelog( "%-35s: %g\n", "UnitVelocity_in_cm_per_s", All.UnitVelocity_in_cm_per_s );
    writelog( "%-35s: %g\n", "Gravity constant", All.G );
    writelog( "%-35s: %g\n", "Hubble", All.Hubble );
    if ( All.MpcFlag != 1 ) {
        All.MpcFlag = 1000;
    }
    g2c.cm       = All.UnitLength_in_cm / header.HubbleParam;
    g2c.g        = All.UnitMass_in_g / header.HubbleParam;
    g2c.s        = All.UnitTime_in_s / header.HubbleParam;
    g2c.erg      = g2c.g * SQR(g2c.cm) / SQR( g2c.s );
    put_block_line;
}

