#include "allvars.h"

#define writelog1( a, b  ) writelog( "%-35s: %g\n", a, b )
void set_units() {

    writelog( "Set Units...\n" );
    All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
    All.UnitDensity_in_cgs = All.UnitMass_in_g / pow( All.UnitLength_in_cm, 3 );
    All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow( All.UnitLength_in_cm,2 ) / pow( All.UnitTime_in_s, 2 );
    All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow( All.UnitTime_in_s, 2 );
    All.G = GRAVITY / pow( All.UnitLength_in_cm, 3 ) * All.UnitMass_in_g * pow( All.UnitTime_in_s, 2 );
    All.Hubble = HUBBLE * All.UnitTime_in_s;

    writelog1( "UnitMass_in_g", All.UnitMass_in_g );
    writelog1( "UnitTime_in_s", All.UnitTime_in_s );
    writelog1( "UnitLength_in_cm", All.UnitLength_in_cm );
    writelog1( "UnitDensity_in_cgs", All.UnitDensity_in_cgs );
    writelog1( "UnitEnergy_in_cgs", All.UnitEnergy_in_cgs );
    writelog1( "UnitPressure_in_cgs", All.UnitPressure_in_cgs );
    writelog1( "UnitVelocity_in_cm_per_s", All.UnitVelocity_in_cm_per_s );
    writelog1( "Gravity constant", All.G );
    writelog1( "Hubble", All.Hubble );

    if ( All.MpcFlag != 1 ) {
        All.MpcFlag = 1000;
    }

    g2c.cm       = All.UnitLength_in_cm / All.HubbleParam;
    g2c.g        = All.UnitMass_in_g / All.HubbleParam;
    g2c.s        = All.UnitTime_in_s / All.HubbleParam;
    g2c.erg      = g2c.g * SQR(g2c.cm) / SQR( g2c.s );

    aux_c.e_mec = ELECTRON_CHARGE / ( ELECTRON_MASS * LIGHT_SPEED );

    put_block_line;
}

