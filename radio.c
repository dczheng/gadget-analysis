#include "allvars.h"

#define TAB_F_N 100000
#define BESSEL_UPPER_LIMIT   ((double)700)
#define F_X_MAX              ((double)600)
#define F_X_MIN              ((double)1e-5)

double *tab_F_V;

double F_integrand( double x, void *params ) {
    if ( x > BESSEL_UPPER_LIMIT )
        return 0;
    return gsl_sf_bessel_Knu( 5.0/3.0, x );
}

void F( double x, double *r, double *err ) {

    gsl_function Func;
    Func.function = &F_integrand;
    Func.params = NULL;

    gsl_integration_qagiu( &Func, x,
            GSL_INTE_ERR_ABS, GSL_INTE_ERR_REL, GSL_INTE_WS_LEN,
            inte_ws, r, err );
    /*
    *err = 1e-1;

    x = 1e-5;
    printf( "x: %g\n", x );

    *r = qtrap( &F_integrand, NULL, x, F_INTE_UPPER_LIMIT, *err );
    */

}

void init_tab_F() {
    double dlogx, x, err, F_v, *buf, err_max, err_max_global;
    int i;

    writelog( "initialize tab_F ...\n" );

    dlogx = log(F_X_MAX/F_X_MIN) / ( TAB_F_N - 1 );

    mymalloc1( tab_F_V, sizeof( double ) * ( TAB_F_N+1 ) );

    err_max = -1;
    for ( i=0; i<=TAB_F_N; i++ ) {
        if ( i % NTask != ThisTask )
            continue;
        x = log( F_X_MIN ) + dlogx * i;
        x = exp( x );
        F( x, &F_v, &err );

        if ( err>err_max )
            err_max = err;

        tab_F_V[i] = log( F_v );
    }

    MPI_Reduce( &err_max, &err_max_global, 1, MPI_DOUBLE,
            MPI_MAX, 0, MPI_COMM_WORLD );

    mymalloc2( buf, sizeof( double ) * ( TAB_F_N+1 ) );
    MPI_Allreduce( tab_F_V, buf, TAB_F_N+1, MPI_DOUBLE,
            MPI_SUM, MPI_COMM_WORLD );
    memcpy( tab_F_V, buf, (TAB_F_N+1) * sizeof( double ) );
    myfree( buf );

    MPI_Barrier( MPI_COMM_WORLD );

    writelog( "initialize tab_F_V ... done.\n" );

}

void free_tab_F() {
    myfree( tab_F_V );
}

double tab_F( double x ) {

    double dx, dlogx, logx, r;
    int i;

    dlogx = log(F_X_MAX/F_X_MIN) / ( TAB_F_N - 1 );
    logx = log( x );
    i = (logx-log(F_X_MIN) ) / dlogx;

    if ( i<0 )
        i = 0;
    if ( i>TAB_F_N-1 )
        i = TAB_F_N-1;

    dx = logx - log(F_X_MIN) - i * dlogx;

    r = exp(tab_F_V[i]) * ( 1-dx ) + exp(tab_F_V[i+1]) * dx;
    //r = exp( tab_F_V[i] );

    return r;

}

struct radio_inte_struct{
    double (*f) ( double, void* );
    void *params;
    double B;
    double nu;
};

double radio_inte( double p, void *params ) {

    double r, x, nu_c;

#ifdef DISABLE_RADIO_F_TAB
        double err;
#endif

    struct radio_inte_struct *ri;

    ri = params;

    nu_c = 0.1875 * ( 1+p*p ) * ri->B * aux_c.e_mec ;  // 0.1875: 3/16

    x = ri->nu / nu_c;

#ifdef DISABLE_RADIO_F_TAB
    F( x, &r, &err );
#else
    r = tab_F( x );
#endif

    r *= x;

    r *= (*(ri->f))( p, ri->params );

    return r;

}

double radio( double (*f)( double, void* ), double *params,
        double B, double nu, double pmin, double pmax ) {

    double r, err, fac, t;
    struct radio_inte_struct ri;

    ri.f = f;
    ri.params = params;
    ri.B = B;
    ri.nu = nu;

    err = 1e-3;

    t = sqrt( nu / BESSEL_UPPER_LIMIT / ( 0.1875 * B * aux_c.e_mec )-1 );

    if ( pmin < t )
        pmin = t;

    r = qtrap( &radio_inte, &ri, pmin, pmax, err );

    fac = sqrt(3) * CUBE( ELECTRON_CHARGE ) * PI / ( 4.0 * ELECTRON_MASS * SQR(LIGHT_SPEED) );
    r *= fac;

    return r;

}

void test_tab_F() {

    double x, dx, xmax, xmin, F_v, tab_F_v, err, err_max, err_mean;
    int i, N;

    xmin = F_X_MIN;
    xmax = F_X_MAX - 1;
    N = 1000;
    dx = log( xmax/xmin ) / ( N-1 );

    /*
    i = 0;
    while( i<N ) {
        x = xmin + i * dx;
        printf( "%g %g\n", x, radio_F_integrand( x, NULL ) );
        i++;
    }
    */

    i = 0;
    err_max = err_mean = 0;

    while( i<N ) {
        x = log(xmin) + i * dx;
        x = exp(x);
        F( x, &F_v, &err );
        tab_F_v = tab_F( x );
        err = fabs( tab_F_v - F_v ) / F_v * 100;
        printf( "x: %g, F: %g, tab_F: %g, err: %.2f%%\n", x, F_v, tab_F_v, err);
        if ( err > err_max )
            err_max = err;
        err_mean += err;
        i++;
        }

        err_mean /= N;

    printf( "err_max: %.2f%%, err_mean: %.2f%%\n", err_max, err_mean );
}

void test_F() {

    double logx_min, logx_max, x, dlogx, F_v, err, err_max, err_mean;
    int logx_N, i;
    FILE *fd;

    logx_min = -5;
    logx_max = 1;
    logx_N = 30;

    dlogx = ( logx_max-logx_min ) / ( logx_N - 1 );

    fd = fopen( "F_x.dat", "w" );

    err_max = err_mean = 0;

    for( i=0; i<logx_N; i++ ) {
        x = logx_min + i * dlogx;
        x = pow( 10, x );
        F( x, &F_v, &err );
        fprintf( fd, "%g %g %g %g\n", x, x*F_v, F_v, err );
        if ( err > err_max )
            err_max = err;

        err_mean += err;
    }

    err_mean /= logx_N;

    printf( "err_max: %g, err_mean: %g\n", err_max, err_mean );

    fclose( fd );

}
