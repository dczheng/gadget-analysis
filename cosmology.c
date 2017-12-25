#include "allvars.h"

double E_a ( double a ) {
    double E2;
    E2 =  header.Omega0 / ( a*a*a ) +
          ( 1-header.Omega0-header.OmegaLambda ) / ( a*a ) +
          header.OmegaLambda;
    return sqrt( E2 );
}

double H_a( double a ) {
    double H;
    H = header.HubbleParam * 100 * E_a( a );
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
    d *= ( LIGHT_SPEED / 1e5 ) / ( header.HubbleParam * 100 ); //Mpc
    return d;
}

double angular_distance( double a ) {
    return comoving_distance( a ) * a;
}

double luminosity_distance( double a ) {
    return comoving_distance( a ) / a;
}
