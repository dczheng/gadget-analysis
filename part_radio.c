#include "allvars.h"

#define TAB_F_N 100000
#define F_INTE_UPPER_LIMIT   ((double)700)
#define F_X_MAX              ((double)600)
#define F_X_MIN              ((double)1e-5)
#define BCMB0 (3.24e-6) // Gauss
//#define BCMB0 (0) // Gauss
//#define DISABLE_TAB_F

double *tab_F_V;


double F_integrand( double x, void *params ) {
    if ( x > F_INTE_UPPER_LIMIT )
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

double radio_integrand( double p, void *params) {

    double *pa, nu_c, B, nu, r, x, a;
#ifdef DISABLE_TAB_F
    double err;
#endif
    pa = params;
    nu = pa[0];
    a = pa[1];
    B = pa[2];
    nu_c = 0.1875 * ( 1+p*p ) * B * aux_c.e_mec ;  // 0.1875: 3/16
    x = nu / nu_c;
#ifdef DISABLE_TAB_F
    F( x, &r, &err );
#else
    r = tab_F( x );
#endif

    r *= x;
    /*
    printf( "nu: %g, nu_c: %g, a: %g, B: %g, p: %g, x: %g, F(x): %g,",
            nu, nu_c, a, B, p, x, r );
            */

    r *= pow( p, -a );

    //printf( "r: %g\n", r );

    return r;

}

double particle_radio2( double nu,  SphParticleData *part ) {

    double r, err, fac, B, c, a, qmin, qmax, params[3];

    //gsl_function Func;

    c = part->CRE_C;
    c = c * part->Density / ( ELECTRON_MASS /(g2c.g) );
    c /= CUBE( g2c.cm );
    a = part->CRE_Alpha;
    qmin = part->CRE_qmin;
    qmax = part->CRE_qmax;


    if ( c * a * qmin * qmax == 0 )
        return 0;

    B = pow( part->B[0], 2 ) + pow( part->B[1], 2 ) + pow( part->B[2], 2 );
    B += SQR( BCMB0 ) * pow( All.Time, -2 );
    B = sqrt( B );

    //Func.function = &radio_integrand;
    //Func.params = params;
    params[0] = nu;
    params[1] = a;
    params[2] = B;


    qmin = nu / F_INTE_UPPER_LIMIT / ( 0.1875 * B * aux_c.e_mec )-1;

    qmin = sqrt( qmin );
    //printf( "qmin: %g\n", qmin );

    fac = c * sqrt(3) * CUBE( ELECTRON_CHARGE ) * PI / ( 4.0 * ELECTRON_MASS * SQR(LIGHT_SPEED) );

    /*
    printf( "cm: %g, g: %g, rho: %g, c: %g, fac: %g\n",
            g2c.cm, g2c.g,
            part->Density, c, fac );
            */

    err = 1e-1;
    r = qtrap( &radio_integrand, params, qmin, qmax, err );

    /*
    size_t neval;
    gsl_integration_qng( &Func, qmin, qmax,
            GSL_INTE_ERR_ABS, 2, &r, &err, &neval );
            */
    /*
    gsl_integration_qag( &Func, qmin, qmax,
            GSL_INTE_ERR_ABS, 1e-5, GSL_INTE_WS_LEN, GSL_INTE_KEY,
            inte_ws, &r, &err );
            */

    //printf( "r: %e\n", r );

    r *= fac;

    r = r * ( 4.0/3.0 * PI * CUBE( All.SofteningTable[0] * g2c.cm ) );

    //printf( "%g\n", r );

    return r;

}

double particle_radio( double nu, long i ) {

    return particle_radio2( nu, &SphP[i] );

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

void test_particle_radio() {

    SphParticleData part;
    double nu_min, nu_max, nu, dlognu, P, qmax_min, qmax_max, dlogqmax;
    int i, N, qmaxn, k;
    FILE *fd;

    part.CRE_C = 1 * CUBE(g2c.cm) * ( ELECTRON_MASS/(g2c.g) );
    part.CRE_Alpha = 2.01;
    part.CRE_qmin = 1;
    part.CRE_qmax = 1e5;
    part.Density = 1;
    part.B[0] = 1e-6;
    part.B[1] = part.B[2] = 0;

    N = 100;
    nu_min = 1e6;
    nu_max = 1e9;
    dlognu = log10( nu_max/nu_min ) / ( N-1 );

    qmaxn = 5;
    qmax_min = 5e3;
    qmax_max = 1e5;
    dlogqmax = log10( qmax_max/qmax_min ) / ( qmaxn-1 );

    fd = fopen( "./particle_radio.dat", "w" );

    fprintf( fd, "0 0 0 0 0 " );

    for( i=0; i<N; i++ ) {
        nu = log10(nu_min) + i * dlognu;
        nu = pow( 10, nu );
        fprintf( fd, "%g ", nu );
    }

    fprintf( fd, "\n" );

    for( k=0; k<qmaxn; k++ ) {
        part.CRE_qmax = log10( qmax_min ) + k*dlogqmax;
        part.CRE_qmax = pow( 10, part.CRE_qmax );

        //part.CRE_qmax = qmax_max;

        fprintf( fd, "%g %g %g %g %g ",
                part.CRE_C /( CUBE(g2c.cm) * ( ELECTRON_MASS/(g2c.g) )),
                part.CRE_Alpha,
                part.CRE_qmin,
                part.CRE_qmax,
                part.B[0] );
        printf(  "%g %g %g %g %g\n",
                part.CRE_C /( CUBE(g2c.cm) * ( ELECTRON_MASS/(g2c.g) )),
                part.CRE_Alpha,
                part.CRE_qmin,
                part.CRE_qmax,
                part.B[0] );

        for( i=0; i<N; i++ ) {
            nu = log10(nu_min) + i * dlognu;
            nu = pow( 10, nu );
            P = particle_radio2( nu, &part );
            //printf( "nu: %g, P: %g\n", nu, P );
            //break;
            fprintf( fd, "%g ", P );
        }

        fprintf( fd, "\n" );

        //break;

    }

    fclose( fd );

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

void test_radio() {


    All.Time = 1;
    All.HubbleParam = 0.7;
    set_units();
    put_block_line;

#ifndef DISABLE_TAB_F
    init_tab_F();
#endif



    if ( ThisTask == 0 ) {
        //test_F();
        test_particle_radio();
    }

    MPI_Barrier( MPI_COMM_WORLD );
#ifndef DISABLE_TAB_F
    free_tab_F();
#endif
    endrun(20181004);

}

/*
void compute_radio( double v ) {

    mymalloc( img, sizeof( double ) * SQR( PicSize ) );
    memset( img, 0, sizeof( double ) * SQR( PicSize ) );

    dy = dx = BoxSize / PicSize;
    h = All.SofteningTable[0] * ( g2c.cm );
    V = 4.0 / 3.0 * PI * pow( h, 3 );
    ang = h / ang_dis / PI * 180.0 * 60;
    beam = pow( 10.0 / 60, 2.0 ) / ( 4.0 * log(2) );
    img_max = DBL_MIN;
    img_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        xi = (int)( P[i].Pos[0] / dx );
        yi = (int)( P[i].Pos[1] / dy );
        tmp = img[xi * PicSize + yi] += SphP[i].P * V / ( 4*PI*pow(lum_dis,2) ) / SQR(ang)
            * beam * 1e25;
        //tmp = img[xi * PicSize + yi] += SphP[i].P * V / ( 4*PI*pow(lum_dis,2) ) / SQR(ang)
        //    * 1e25;
        img_max = ( tmp > img_max ) ? tmp : img_max;
        img_min = ( tmp < img_min && tmp > 0 ) ? tmp : img_min;
    }
}
*/
/*
double particle_radio( double nu, long i ) {

    C = SphP[i].CRE_C;
    //printf( "%g\n", C );
    C = C * SphP[i].Density / ( ELECTRON_MASS /(g2c.g) );
    C /= CUBE( g2c.cm );

    B = sqrt( pow( SphP[i].B[0], 2 ) + pow( SphP[i].B[1], 2 ) + pow( SphP[i].B[2], 2 ) );
    Ub = B * B / ( 8 * PI );
    Ub += SQR( BCMB0 ) * pow( All.Time, -4 ) / ( 8*PI );

    nuL = ELECTRON_CHARGE * B / ( 2 * PI * ELECTRON_MASS * LIGHT_SPEED );

    if ( sqrt( nu/nuL-1 ) < SphP[i].CRE_qmin ||
         sqrt( nu/nuL-1 ) > SphP[i].CRE_qmax )
        return 0;

    P = C * 2.0 / 3.0 * LIGHT_SPEED * Ub * THOMSON_CROSS_SECTION *
        pow( nu / nuL-1, (1-SphP[i].CRE_Alpha) / 2 ) / nuL;

    //printf( "%g\n", P );
    //
    P = P * ( 4.0/3.0 * PI * CUBE( All.SofteningTable[0] * g2c.cm ) );

    // Unit: erg / Hz / s
    //
    //printf( "Ub: %g, nuL: %g, P: %g\n",  Ub, nuL, P );

    printf( "P: %g, P2: %g\n", P, P2 );

    endrun( 20181005 );

    return P;

}
*/

