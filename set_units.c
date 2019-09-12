#include "allvars.h"

#define writelog1( a ) writelog( "%-35s: %g\n", #a, a )
#define writelog2( a ) writelog( "%-6s: %g\n", #a, a )

void set_units() {

    writelog( "Set Units...\n" );
    UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
    UnitDensity_in_cgs = All.UnitMass_in_g / pow( All.UnitLength_in_cm, 3 );
    UnitEnergy_in_cgs = All.UnitMass_in_g * pow( All.UnitLength_in_cm,2 ) / pow( UnitTime_in_s, 2 );
    UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow( UnitTime_in_s, 2 );
    Hubble = HUBBLE * UnitTime_in_s;
    G = GRAVITY / pow( All.UnitLength_in_cm, 3 ) * All.UnitMass_in_g * pow( UnitTime_in_s, 2 );

    writelog1( All.UnitMass_in_g );
    writelog1( All.UnitLength_in_cm );
    writelog1( All.UnitVelocity_in_cm_per_s );
    writelog1( UnitTime_in_s );
    writelog1( UnitDensity_in_cgs );
    writelog1( UnitEnergy_in_cgs );
    writelog1( UnitPressure_in_cgs );
    writelog1( Hubble );
    writelog1( G );

    if ( All.MpcFlag != 1 ) {
        All.MpcFlag = 1000;
    }

    g2c.cm       = All.UnitLength_in_cm / HubbleParam;
    g2c.g        = All.UnitMass_in_g    / HubbleParam;
    g2c.s        = UnitTime_in_s        / HubbleParam;
    g2c.erg      = g2c.g * SQR(g2c.cm)  / SQR(g2c.s);

    cuc.m_e = ELECTRON_MASS;
    cuc.m_p = PROTONMASS;
    cuc.e = ELECTRON_CHARGE;
    cuc.c = LIGHT_SPEED;
    cuc.sigma_t = THOMSON_CROSS_SECTION;
    cuc.e_mec = cuc.e / ( cuc.m_e * cuc.c );
    cuc.G = GRAVITY;
    cuc.c2 = cuc.c * cuc.c;
    cuc.mec2 = cuc.m_e * cuc.c2;

    /*
    writelog( "cgs: \n" );
    writelog2( cuc.m_e );
    writelog2( cuc.m_p );
    writelog2( cuc.e );
    writelog2( cuc.c );
    writelog2( cuc.e_mec );
    writelog2( cuc.G );
    writelog2( cuc.c2 );
    writelog2( cuc.mec2 );
    */

    guc.m_e = cuc.m_e / g2c.g;
    guc.m_p = cuc.m_p / g2c.g;
    guc.c = cuc.c / ( g2c.cm / g2c.s );
    guc.e = cuc.e / ( pow( g2c.cm, 1.5 ) * sqrt(g2c.g) / g2c.s );
    guc.e_mec = guc.e / ( guc.m_e * guc.c );
    guc.c2 = guc.c * guc.c;
    guc.mec2 = guc.m_e * guc.c2;
    guc.G = G;
    guc.sigma_t = cuc.sigma_t / ( SQR(g2c.cm) );

    SofteningTable[0] = All.SofteningGas;
    SofteningTable[1] = All.SofteningHalo;
    SofteningTable[2] = All.SofteningDisk;
    SofteningTable[3] = All.SofteningBulge;
    SofteningTable[4] = All.SofteningStar;
    SofteningTable[5] = All.SofteningBndry;
    Start[0] = All.StartX;
    Start[1] = All.StartY;
    Start[2] = All.StartZ;
    End[0] = All.EndX;
    End[1] = All.EndY;
    End[2] = All.EndZ;
    PicSize = All.PicSize;

    /*
    writelog( "gadget: \n" );
    writelog2( guc.m_e );
    writelog2( guc.m_p );
    writelog2( guc.e );
    writelog2( guc.c );
    writelog2( guc.e_mec );
    writelog2( guc.G );
    writelog2( guc.c2 );
    writelog2( guc.mec2 );
    */
}

