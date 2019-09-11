#include "stdio.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf_bessel.h"
//#include "../mymath.c"

#define GSL_BESSEL_NU (5.0/3.0)
#define GSL_BESSEL_UPPER_LIMIT (700)
#define GSL_BESSEL_BREAKPOINT (10)

double my_gsl_bessel_Knu( double x, void *p ) {

    if ( x > GSL_BESSEL_UPPER_LIMIT )
        return 0;

     return gsl_sf_bessel_Knu( GSL_BESSEL_NU, x );

}

void gsl_bessel_Knu() {

    double xmin, xmax, dlogx, x, r, rr, err, rrr;
    int i, N;
    FILE *fd;
    gsl_function F;
    gsl_integration_workspace *inte_ws;

    xmin = 1e-8;
    xmax = 600;
    N = 2000;

    F.function = &my_gsl_bessel_Knu;
    F.params = NULL;
    inte_ws = gsl_integration_workspace_alloc( 10000 );

    dlogx = log( xmax / xmin ) / ( N-1 );

    fd = fopen( "gsl_bessel_Knu.dat", "w" );

    gsl_integration_qagiu( &F, GSL_BESSEL_BREAKPOINT,
                0, 1e-2, 10000, inte_ws, &rrr, &err );

    for( i=0; i<N; i++ ) {

        x = exp( log( xmin ) + i * dlogx );
        r = my_gsl_bessel_Knu( x, NULL );

        if ( x < GSL_BESSEL_BREAKPOINT ) {
            gsl_integration_qag( &F, x, GSL_BESSEL_BREAKPOINT,
                0, 1e-2, 10000,
                GSL_INTEG_GAUSS61,
                inte_ws, &rr, &err );

            rr = rr + rrr;
        }
        else
            gsl_integration_qagiu( &F, x,
                0, 1e-2, 10000, inte_ws, &rr, &err );

        /*
        rrr = qtrap( &my_gsl_bessel_Knu, NULL, x, GSL_BESSEL_UPPER_LIMIT, 1e-2 );
        */

        fprintf( fd, "%g %g %g\n", x, r,rr );

    }

    gsl_integration_workspace_free( inte_ws );

    fclose( fd );

}

int main() {

    gsl_bessel_Knu();

    system( "./plot_gsl_test.py" );
    return 0;

}
