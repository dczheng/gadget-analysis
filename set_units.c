#include "allvars.h"

#define writelog1( a ) writelog( "%-35s: %g\n", #a, a )
void set_units() {

    put_header( "Set Units" );
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
    writelog1( HubbleParam );

    if ( All.MpcFlag != 1 ) {
        All.MpcFlag = 1000;
    }

    g2c.cm       = All.UnitLength_in_cm / HubbleParam;
    g2c.g        = All.UnitMass_in_g    / HubbleParam;
    g2c.s        = UnitTime_in_s        / HubbleParam;

    g2c.density  = g2c.g / CUBE( g2c.cm );
    g2c.erg      = g2c.g * SQR(g2c.cm)  / SQR(g2c.s);
    g2c.pressure = g2c.g / ( g2c.cm * SQR(g2c.s) );

    cuc.m_e = ELECTRON_MASS;
    cuc.m_p = PROTONMASS;
    cuc.e = ELECTRON_CHARGE;
    cuc.c = LIGHT_SPEED;
    cuc.sigma_t = THOMSON_CROSS_SECTION;
    cuc.e_mec = cuc.e / ( cuc.m_e * cuc.c );
    cuc.G = GRAVITY;
    cuc.c2 = cuc.c * cuc.c;
    cuc.mec2 = cuc.m_e * cuc.c2;

    writelog( "cgs: \n" );
    writelog1( cuc.m_e );
    writelog1( cuc.m_p );
    writelog1( cuc.e );
    writelog1( cuc.c );
    writelog1( cuc.sigma_t );
    writelog1( cuc.e_mec );
    writelog1( cuc.G );
    writelog1( cuc.c2 );
    writelog1( cuc.mec2 );

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
    PicSize = All.PicSize;

    writelog( "gadget: \n" );
    writelog1( guc.m_e );
    writelog1( guc.m_p );
    writelog1( guc.e );
    writelog1( guc.c );
    writelog1( cuc.sigma_t );
    writelog1( guc.e_mec );
    writelog1( guc.G );
    writelog1( guc.c2 );
    writelog1( guc.mec2 );


    All.FreqMin *= 1e6 / Time;
    All.FreqMax *= 1e6 / Time;
    All.GroupRadFreq *= 1e6 / Time;
    All.GroupRadFreq1 *= 1e6 / Time;
    All.OutputGroupFreq *= 1e6 / Time;
    All.RadSliceFreq *= 1e6 / Time;
    All.GroupRadProfileFreq *= 1e6 / Time;

    writelog1( All.FreqMin );
    writelog1( All.FreqMax );
    writelog1( All.GroupRadFreq );
    writelog1( All.GroupRadFreq1 );
    writelog1( All.GroupRadFreq );
    writelog1( All.RadSliceFreq );
    writelog1( All.GroupRadProfileFreq );

}

