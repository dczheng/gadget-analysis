#include "allvars.h"

double trapzd( double (*func)( double, void* ),
        void *params, double a, double b, int n ) {

    int i, num;
    double sum, step, x;
    static double s;

    if ( n == 1 ) {
        s = 0.5 * ( (*func)( a, params ) + (*func)(b, params) ) * (b-a);
        return s;
    }

    for( num=1, i=1; i<n-1; i++ ) num <<= 1;

    step = (b-a) / num;

    x = a+0.5*step;

    for( sum=0, i=0; i<num; i++, x+=step )
        sum += (*func)(x, params);

    s = 0.5 * ( s + (b-a)*sum / num );

    return s;

}

#define REFINE_MAX 20
double qtrap( double (*func)( double, void* ),
        void *params, double a, double b, double err ) {

    double s, olds;
    int i;

    olds = -DBL_MAX;

    for( i=1; i<=REFINE_MAX; i++ ) {

        s = trapzd( func, params, a, b, i );

        if ( i>5 )
            if ( fabs( s-olds ) < err*fabs(olds) ||
                    ( s == 0 && olds == 0 ) )
                return s;

        olds = s;
    }

    printf( "a: %g, b: %g\n", a, b );
    printf( "Too many steps in qtrap!\n" );
    RAISE_SIGSTOP();
    endrun( 20181008 );
     return  s;

}
