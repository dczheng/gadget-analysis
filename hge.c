#include "allvars.h"

#define F( c, a, q1, q2, p )                 ( ( p<q1 || p>q2 ) ? 0 : c * pow( p, -a ) )

#define BETA(a, q)                           beta( (a-2)*0.5, (3-a)*0.5, 1/(1+q*q) )

#define NUMBER_DENSITY( c, a, q )            ( c * pow( q, 1.0-a ) / ( a-1.0 ) )

#define ENERGY( c, a, q )                    ( LS2 * c / (a-1.0) * ( 0.5 * BETA(a, q) + (sqrt( 1+q*q ) - 1.0 ) * pow( q, 1-a ) ) )
#define ENERGY2( c, a, q1, q2 )              ( ENERGY0( c, a, q1 ) - ENERGY0( c, a, q2 ) )

#define MEAN_ENERGY( a, q )                  ( ( 0.5 * pow( q, a-1.0 ) * BETA(a, q) + sqrt( 1+q*q ) - 1.0 ) * MEC2 )
#define MEAN_ENERGY2( a, q1, q2 )            ( ENERGY02( 1, a, q1, q2 ) / NUMBER_DENSITY2( 1, a, q1, q2 ) * ELECTRON_MASS )

#define PRESSURE( c, a, q )                  ( LS2 * c / 6 * BETA( a, q ) )
#define PRESSURE2( c, a, q1, q2 )            ( PRESSURE( c, a, q1 ) - PRESSURE( c, a, q2 ) )


double beta_inte( double x, void *params ) {

    double *p = params;
    return pow( x, p[0]-1.0 ) * pow( 1.0-x, p[1]-1.0 );

}

double beta( double a, double b, double x ) {

    double r, err;
    gsl_function F;
    double params[2];

    if ( a<0 ) {
        RAISE_SIGSTOP();
        endrun( 20190922 );
    }

    F.function = &beta_inte;

    if ( b>0.0 )
        return gsl_sf_beta_inc( a, b, x ) * gsl_sf_beta( a, b );

    F.params = params;
    params[0] = a;
    params[1] = b;

    gsl_integration_qag( &F, 0, x,
            GSL_INTE_ERR_ABS, GSL_INTE_ERR_REL, GSL_INTE_WS_LEN, GSL_INTE_KEY,
            inte_ws, &r, &err );

    return r;

}

void hge_pressure_pdf() {

}
