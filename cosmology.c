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

void compute_cosmo_quantities() {
    writelog( "compute cosmology quantities ...\n" );
    All.Time = 1 / ( All.RedShift + 1 );
    All.Hubble_a = hubble_function( All.Time );
    All.RhoCrit = 3 * SQR( All.Hubble_a ) / ( 8*PI*All.G );
    All.RhoBaryon = All.OmegaBaryon * All.RhoCrit;
    writelog( "%-35s: %g\n", "Time", All.Time );
    writelog( "%-35s: %g\n", "Hubble_a", All.Hubble_a );
    writelog( "%-35s: %g\n", "RhoBaryon", All.RhoBaryon );
    writelog( "%-35s: %g\n", "RhoCrit", All.RhoCrit );
    writelog( "compute cosmology quantities ... done.\n" );
    put_block_line;
}

