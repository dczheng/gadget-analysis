#include "allvars.h"

#define BETA(a, q)                           beta( (a-2)*0.5, (3-a)*0.5, 1/(1+q*q) )

#define NUMBER_DENSITY( c, a, q )            ( (c) * pow( (q), 1.0-(a) ) / ( (a)-1.0 ) )
#define NUMBER_DENSITY2( c, a, q1, q2 )      ( NUMBER_DENSITY( (c), (a), (q1) ) - NUMBER_DENSITY( (c), (a), (q2) ) )

#define ENERGY( c, a, q )                    ( guc.c2 * (c) / ((a)-1.0) * ( 0.5 *  BETA((a), (q)) + (sqrt( 1+(q)*(q)) - 1.0 ) * pow( (q), 1-(a) ) ) )
#define ENERGY2( c, a, q1, q2 )              ( ENERGY( (c), (a), (q1) ) -  ENERGY( (c), (a), (q2) ) )

#define MEAN_ENERGY( a, q )                  ( ( 0.5 * pow( (q), (a)-1.0 ) * BETA((a), (q)) + sqrt( 1+((q)*(q)) ) - 1.0 ) * guc.mec2 )
#define MEAN_ENERGY2( a, q1, q2 )            ( ENERGY2( 1, (a), (q1), (q2) ) / NUMBER_DENSITY2( 1, (a), (q1), (q2) ) * guc.m_e )


#define PRESSURE( c, a, q )                  ( guc.c2 * (c) / 6 * BETA( (a), (q) ) )
#define PRESSURE2( c, a, q1, q2 )            ( PRESSURE( (c), (a), (q1) ) - PRESSURE( (c), (a), (q2) ) )

#define cre_pressure( i )                    ( PRESSURE2( SphP[i].CRE_C, SphP[i].CRE_Alpha, SphP[i].CRE_qmin, SphP[i].CRE_qmax ) )


double beta_inte( double x, void *params ) {

    double *p = params;
    return pow( x, p[0]-1.0 ) * pow( 1.0-x, p[1]-1.0 );

}

double beta( double a, double b, double x ) {

    double r, err;
    gsl_function F;
    double params[2];

    /*
    gsl_integration_workspace *inte_ws_local;
    inte_ws_local = gsl_integration_workspace_alloc( GSL_INTE_WS_LEN );
    */

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

    //gsl_integration_workspace_free( inte_ws_local );

    return r;

}

void compute_cre_pressure() {

#if defined(CREPPDF) || defined(CRPPDF)

    long i, t;
    t = N_Gas / 10;
    mytimer_start();
    writelog( "compute cre pressure ...\n" );

    for( i=0; i<N_Gas; i++ ) {
        if ( i % t == 0 )
            writelog( "[%10li] [%10li] [%5.1f%%]\n", i, N_Gas, ( (double)(i) ) / N_Gas * 100 );
        if ( i % NTask_Local != ThisTask_Local )
            continue;
        if ( SphP[i].CRE_n != 0 ) {
            SphP[i].CRE_P = cre_pressure( i ) * (SphP[i].Density / CUBE(Time));
        }
    }
    mytimer_end();
#endif

}

