#include "allvars.h"

#define TAB_RADIO_F_N 1000000
#define RADIO_F_INTE_UPPER_LIMIT   ((double)700)
#define RADIO_F_X_MAX              ((double)600)
#define RADIO_F_X_MIN              ((double)1e-6)
#define BCMB0 (3.24e-6) // Gauss

double *tab_radio_F_F;


double radio_F_integrand( double x, void *params ) {
    return gsl_sf_bessel_Knu( 5.0/3.0, x );
}

void radio_F( double x, double *r, double *err ) {

    gsl_function F;
    F.function = &radio_F_integrand;
    F.params = NULL;

    if ( x >= RADIO_F_INTE_UPPER_LIMIT ) {
        *r = 0;
        *err = 0;
        return;
    }

    gsl_integration_qag( &F, x, RADIO_F_INTE_UPPER_LIMIT,
            GSL_INTE_ERR_ABS, GSL_INTE_ERR_REL, GSL_INTE_WS_LEN, GSL_INTE_KEY,
            inte_ws, r, err );

}

void init_tab_radio_F() {
    double dlogx, x, err, F, *buf;
    int i;
    FILE *fd;
    char *fn_tab = "./tab_F.dat";

    writelog( "initialize tab_radio_F_F ...\n" );
    dlogx = log(RADIO_F_X_MAX/RADIO_F_X_MIN) / ( TAB_RADIO_F_N - 1 );

    mymalloc1( tab_radio_F_F, sizeof( double ) * ( TAB_RADIO_F_N+1 ) );

    if ( access( fn_tab, 0 ) != -1 ) {
        writelog( "read tab_F.dat ...\n" );

        if ( ThisTask == 0 ) {
            fd = fopen( fn_tab, "r" );
            for ( i=0; i<=TAB_RADIO_F_N; i++ )
                fscanf( fd, "%lf\n", &tab_radio_F_F[i] );
            fclose( fd );
        }

        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Bcast( tab_radio_F_F, (TAB_RADIO_F_N+1) * sizeof( double ),
                MPI_BYTE, 0, MPI_COMM_WORLD );
        return;
    }

    for ( i=0; i<=TAB_RADIO_F_N; i++ ) {
        if ( i % NTask != ThisTask )
            continue;
        x = log( RADIO_F_X_MIN ) + dlogx * i;
        x = exp( x );
        radio_F( x, &F, &err );
        tab_radio_F_F[i] = log( F );
    }

    mymalloc2( buf, sizeof( double ) * ( TAB_RADIO_F_N+1 ) );
    MPI_Allreduce( tab_radio_F_F, buf, TAB_RADIO_F_N+1, MPI_DOUBLE,
            MPI_SUM, MPI_COMM_WORLD );
    memcpy( tab_radio_F_F, buf, (TAB_RADIO_F_N+1) * sizeof( double ) );
    myfree( buf );

    writelog( "write tab_F.dat ...\n" );

    if ( ThisTask == 0 ) {
        fd = fopen( fn_tab, "w" );
        for ( i=0; i<=TAB_RADIO_F_N; i++ )
            fprintf( fd, "%lf\n", tab_radio_F_F[i] );
        fclose( fd );
    }

    MPI_Barrier( MPI_COMM_WORLD );

    writelog( "initialize tab_radio_F_F ... done.\n" );

}

void free_tab_radio_F() {
    myfree( tab_radio_F_F );
}

double tab_radio_F( double x ) {

    double dx, dlogx, logx, r;
    int i;

    dlogx = log(RADIO_F_X_MAX/RADIO_F_X_MIN) / ( TAB_RADIO_F_N - 1 );
    logx = log( x );
    i = (logx-log(RADIO_F_X_MIN) ) / dlogx;

    if ( i<0 )
        i = 0;
    if ( i>TAB_RADIO_F_N-1 )
        i = TAB_RADIO_F_N-1;

    dx = logx - log(RADIO_F_X_MIN) - i * dlogx;

    r = exp(tab_radio_F_F[i]) * ( 1-dx ) + exp(tab_radio_F_F[i+1]) * dx;
    //r = exp( tab_radio_F_F[i] );

    return r;

}

#define DIS_FUNC ( C, a, pmin, pmax, p ) ( ( p<pmin || p>pmax ) ? 0 : C * pow( p, -a ) )

double radio_integrand( double p, void *params) {

    double *pa, nu_c, B, nu, r, x, a;
    pa = params;
    nu = pa[0];
    a = pa[1];
    B = pa[2];
    nu_c = 0.1875 * ( 1+p*p ) * B * aux_c.e_mec ;  // 0.1875: 3/16
    x = nu / nu_c;
    r = tab_radio_F( x ) * x;
    r *= pow( p, -a );
    //printf( "%g\n", r );

    return r;

}

double particle_radio2( double nu,  SphParticleData *part ) {

    double r, err, fac, B, c, a, qmin, qmax, params[3];

    gsl_function F;

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

    /*
    nu = 1e+8;
    qmin = 10;
    B = 1e-5;
    a = 5;
    qmax = 96.961304;
    printf( "nu: %g, c: %g, a: %10f, q: ( %10f, %10f ), B: %10f\n", nu, c, a, qmin, qmax, B );
    */

    F.function = &radio_integrand;
    F.params = params;
    params[0] = nu;
    params[1] = a;
    params[2] = B;


    fac = c * sqrt(3) * CUBE( ELECTRON_CHARGE ) * PI / ( 4.0 * ELECTRON_MASS * SQR(LIGHT_SPEED) );

    //printf( "fac: %g\n", fac );

    gsl_integration_qag( &F, qmin, qmax,
            GSL_INTE_ERR_ABS, 1e-2, GSL_INTE_WS_LEN, GSL_INTE_KEY,
            inte_ws, &r, &err );

    r *= fac;

    r = r * ( 4.0/3.0 * PI * CUBE( All.SofteningTable[0] * g2c.cm ) );

    //printf( "%g\n", r );

    return r;

}

double particle_radio( double nu, long i ) {

    double C, B, Ub, nuL, P, P2;

    P2 = particle_radio2( nu, &SphP[i] );

    return P2;

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

void test_F() {

    double x, dx, xmax, xmin, F, tab_F, err, err_max, err_mean;
    int i, N;

    xmin = RADIO_F_X_MIN;
    xmax = RADIO_F_X_MAX - 1;
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
        radio_F( x, &F, &err );
        tab_F = tab_radio_F( x );
        err = fabs( tab_F - F ) / F * 100;
        printf( "x: %g, F: %g, tab_F: %g, err: %.2f%%\n", x, F, tab_F, err);
        if ( err > err_max )
            err_max = err;
        err_mean += err;
        i++;
        }

        err_mean /= N;

    printf( "err_max: %.2f%%, err_mean: %.2f%%\n", err_max, err_mean );
}

void test_calc_radio() {

    SphParticleData part;
    double nu_min, nu_max, nu, dlognu, P;
    int i, N;
    FILE *fd;

    All.Time = 1;
    part.CRE_C = 1;
    part.CRE_Alpha = 2.8;
    part.CRE_qmin = 1;
    part.CRE_qmax = 1e5;
    part.B[0] = 1e-6;
    part.B[1] = part.B[2] = 0;
    N = 100;
    nu_min = 1e7;
    nu_max = 1e9;

    dlognu = log10( nu_max/nu_min ) / ( N-1 );

    fd = fopen( "./test_calc_radio.dat", "w" );

    for( i=0; i<N; i++ ) {
        nu = log10(nu_min) + i * dlognu;
        nu = pow( 10, nu );
        P = particle_radio2( nu, &part );
        fprintf( fd, "%g %g\n", nu, P );
        printf( "%g %g\n", nu, P );
    }

    fclose( fd );

}

void test_radio() {

    init_tab_radio_F();

    if ( ThisTask == 0 ) {
        test_F();
        //test_calc_radio();
    }

    MPI_Barrier( MPI_COMM_WORLD );
    free_tab_radio_F();
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

