#include "allvars.h"

double trapzd_log( double (*func)( double, void* ),
        void *params, double a, double b, int n ) {

    int i, num;
    double sum, step, x;
    static double s;


    if ( n == 1 ) {
        s = 0.5 * ( (*func)( a, params ) + (*func)(b, params) ) * (b-a);
        return s;
    }

    a = log(a);
    b = log(b);

    for( num=1, i=1; i<n-1; i++ ) num <<= 1;

    step = (b-a) / num;

    x = a+0.5*step;

    for( sum=0, i=0; i<num; i++, x+=step )
        sum += (*func)(exp(x), params) * exp(x);

    s = 0.5 * ( s + (b-a)*sum / num );

    return s;

}

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

#define REFINE_MAX 30
double qtrap( double (*func)( double, void* ),
        void *params, int logflag, double a, double b, double err ) {

    double s, olds;
    int i;
    double (*ftrapzd)( double (*)(double, void*),
                void*, double,  double, int );

    if (logflag) {
        if ( a<0 || b<0 ) {
            endruns( "a<0 || b<0 for log scale!\n" );
        }
        ftrapzd = &trapzd_log;
    }
    else {
        ftrapzd = &trapzd;
    }

//#define QTRAP_DEBUB
#ifdef QTRAP_DEBUB
printf( "a: %g, b: %g, err: %g\n", a, b, err );
#endif

    olds = -DBL_MAX;

    for( i=1; i<=REFINE_MAX; i++ ) {

        s = (*ftrapzd)( func, params, a, b, i );

        if ( i>5 )
            if ( fabs( s-olds ) < err*fabs(olds) ||
                    ( s == 0 && olds == 0 ) )
                return s;

        olds = s;
#ifdef QTRAP_DEBUB
        printf( "[%2i] s: %g\n", i, s );
#endif
    }

    printf( "a: %g, b: %g\n", a, b );
    printf( "Too many steps in qtrap!\n" );
    RAISE_SIGSTOP();
    endrun( 20181008 );

    return  s;

}

//#define TEST_MATH

#ifdef TEST_MATH

double f( double x, void *params) {
    //return x*x;
    //return exp(-x*x*1000);
    return exp(-100*x);
}
#endif

void test_math() {
#ifdef TEST_MATH
    double a, b, err, r1, r2;
    a = 1e-5;
    b = 10;
    err = 1e-6;
    if ( ThisTask == 0 ) {
        r1 = qtrap( &f, NULL, 0, a, b, err );
        r2 = qtrap( &f, NULL, 1, a, b, err );
        printf( "no-log: %g, log: %g\n", r1, r2 );
    }
    do_sync( "" );
    endruns( "test-math" );
#else
    return;
#endif
}

double get_solid_angle( double a, double b, double d ) {
    return 4 * asin( a*b / sqrt( (a*a+4*d*d)*(b*b+4*d*d) ) );
}

double get_solid_angle2( double a, double b, double d ) {
    double alpha, beta;
    alpha = a / ( 2 * d );
    beta = b / ( 2 * d );
    return 4 * acos(
                    sqrt(
                            (1+SQR(alpha)+SQR(beta))
                            / ((1+SQR(alpha)) * (1+SQR(beta)))
                        )
                    );
}

