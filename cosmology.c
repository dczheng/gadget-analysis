#include "allvars.h"

double E_a ( double a ) {
    double E2;
    E2 =  All.Omega0 / ( a*a*a ) +
          ( 1-All.Omega0-header.OmegaLambda ) / ( a*a ) +
          All.OmegaLambda;
    return sqrt( E2 );
}

double hubble_function( double a ) {
    double H;
    H = All.Hubble * E_a( a );
    return H;
}

double com_integ( double a, void *params ) {
    return  1 / ( a * a * E_a( a ) );
}

double comoving_distance( double a ) {
    double d, err;
    gsl_function F;
    F.function = &com_integ;
    gsl_integration_qag( &F, a, 1,
            GSL_INTE_ERR_ABS, GSL_INTE_ERR_REL, GSL_INTE_WS_LEN,
            GSL_INTE_KEY, inte_ws, &d, &err );
    d *= ( LIGHT_SPEED / 1e5 ) / ( All.Hubble );
    return d;
}

double angular_distance( double a ) {
    return comoving_distance( a ) * a;
}

double luminosity_distance( double a ) {
    return comoving_distance( a ) / a;
}


#define writelog1( a, b ) writelog( "%-35s: %g\n", a, b )
void compute_cosmo_quantities() {
    writelog( "compute cosmology quantities ...\n" );

    All.Time = header.time;
    All.Time2 = SQR( All.Time );
    All.Time3 = CUBE( All.Time );
    All.Hubble_a = hubble_function( All.Time );
    All.RhoCrit = 3 * SQR( All.Hubble_a ) / ( 8*PI*All.G );
    All.RhoBaryon = All.OmegaBaryon * All.RhoCrit;

    if ( All.Time > 1e-10 ){
        All.ComDis = comoving_distance( All.Time );
        All.AngDis = angular_distance( All.Time );
        All.LumDis = luminosity_distance( All.Time );
    }
    else
        All.ComDis = All.AngDis = All.LumDis = 0;

    writelog1( "Time", All.Time );
    writelog1( "Time2", All.Time2 );
    writelog1( "Time3", All.Time3 );
    writelog1( "Comoving Distance", All.ComDis );
    writelog1( "Angular Distance", All.AngDis );
    writelog1( "Luminosity Distance", All.LumDis );
    writelog1( "Hubble_a", All.Hubble_a );
    writelog1( "RhoBaryon", All.RhoBaryon );
    writelog1( "RhoCrit", All.RhoCrit );

    writelog( "compute cosmology quantities ... done.\n" );
    put_block_line;
}

