#include "allvars.h"

sigma_struct global_ss;

double E_a ( double a ) {

    double E2;
    E2 =  All.Omega0 / ( a*a*a ) +
          ( 1-All.Omega0-header.OmegaLambda ) / ( a*a ) +
          All.OmegaLambda;
    return sqrt( E2 );

}

double hubble_function( double a ) {

    return All.Hubble * E_a( a );

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

double OmegaM( double a ) {

    // only for Omega_r = Omega_k = 0
    // C.P. Dullemond
    return  All.Omega0 /
        ( a + All.Omega0 * ( 1-a )  + All.OmegaLambda * ( CUBE(a)-a ) );

}

double OmegaLambda( double a ) {
    // only for Omega_r = Omega_k = 0
    // C.P. Dullemond
    return  All.OmegaLambda * CUBE(a) /
        ( a + All.Omega0 * ( 1-a )  + All.OmegaLambda * ( CUBE(a)-a ) );

}

double growth_factor( double a ){
    // C.P. Dullemond
    //
    double Om, OL;
    Om = OmegaM( a );
    OL = OmegaLambda( a );
    return 2.5 * Om / ( pow( Om, 4/7 ) - OL + ( 1+0.5*Om ) * ( 1+1/70.0*OL ) );

}

double growth_factor0( double a ) {

    return growth_factor( a ) / growth_factor( 1 );

}

double PowerSpec_Efstathiou(double k) {

  double AA, BB, CC, nu, ShapeGamma;

  ShapeGamma = 0.21;
  AA = 6.4 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  BB = 3.0 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  CC = 1.7 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  nu = 1.13;

  return k / pow(1 + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);

}

double top_hat_filter( double k, void *params ) {

    double kR, w, *R;

    R = params;

    kR = *R * k;


    if ( kR < 1e-4 )
        w = 1;
    else
        w = 3.0 * ( sin(kR) / CUBE( kR ) - cos( kR ) / SQR(kR) );

    //printf(  "R: %g, k: %g, kR: %g, w: %g\n", *R, k, kR, w );
    return w;

}

void top_hat_filter_k_limit( double *k0, double *k1, void *params ) {

    double *R;

    R = params;

    *k0 = 1e-99 / *R;
    *k1 = 400.0 / *R;

}

double dsigma_dk( double k, void *params ) {

    double p, kR, w;

    sigma_struct *ss;
    ss = params;

    p = (*(ss->P))( k );

    w =(*(ss->filter))( k, ss->filter_params );

    //printf( "k: %e, p: %e, w: %e\n", k, p, w );

    return SQR(k) * p * SQR(w);

}

double sigma( sigma_struct ss ) {

    double r, err, k0, k1;

    gsl_function F;
    F.function = &dsigma_dk;
    F.params = &ss;
    ss.filter_k_limit( &k0, &k1, ss.filter_params );

    /*
    for ( k0=1e-10; k0<1e2; k0*=5 )
        printf( "k: %e P: %e\n", k0, PowerSpec_Efstathiou( k0 ) );
    printf( "k0: %g, k1: %g, dsigma_dk: %g, %g\n", k0, k1,
            dsigma_dk( k0, &ss ), dsigma_dk( k1, &ss ) );
            */

    //r = qtrap( &dsigma_dk, &ss, k0, k1, 0.01 );
    gsl_integration_qag( &F, k0, k1,
            GSL_INTE_ERR_ABS, 1e-3, GSL_INTE_WS_LEN,
            GSL_INTE_KEY, inte_ws, &r, &err );
    //printf( "err: %g\n", err );

    return sqrt( r / ( 2 * SQR(M_PI) ) );

}

void test_sigma() {

    sigma_struct ss;
    double R;

    printf( "test sigma ...\n" );

    R = 8000.0;
    ss.filter = &top_hat_filter;
    ss.filter_params = &R;
    ss.filter_k_limit = &top_hat_filter_k_limit;
    ss.P = &PowerSpec_Efstathiou;

    printf( "Sigma8: %g\n", sigma( ss ) );

    endrun( 20181025 );

}

double top_hat_dsigma2dmdk( double k, void *params) {

    double p, w, dwdr, kR, *R, kR3, kR2, drdm, R2, d2fac;

    sigma_struct *ss;
    ss = params;

    R = ss->filter_params;

    p = (*(ss->P))( k );

    R2 = *R * *R;
    kR = *R * k;
    kR2 = kR * kR;
    kR3 = kR2 * kR;

    w = (*(ss->filter))( k, ss->filter_params );

    if ( kR < 1e-10 )
        dwdr = 0;
    else
        dwdr = 9 * cos( kR ) * k / kR3 + 3.0 * sin(kR) * ( 1-3/kR2 ) / (kR* *R);

    drdm = 1 / ( 4 * M_PI * All.RhoM * R2 );

    return 0;
    //return SQR(k) * p * 2 * w * dwdr * drdm * d2fac;

}

double dsigma2dm( double M, sigma_struct ss ) {

    /*
    double R, d2fac, k0, k1, r, err, sigmaM;

    R = 3 * M / pow( 3 * M_PI * All.RhoM, 1.0/3.0 );

    ss.filter = &top_hat_filter;
    ss.filter_params = &R;
    ss.filter_k_limit = &top_hat_filter_k_limit;
    ss.norm = All.SigmaNorm;
    sigmaM = sigma( ss );

    ss.dsigma2dmdk = &top_hat_dsigma2dmdk;
    d2fac = M*1e4 / sigmaM;
    ss.dsigma2dmdk_params = &d2fac;

    gsl_function F;
    F.function = ss.dsigma2dmdk;
    F.params = ss.dsigma2dmdk_params;

    ss.filter_k_limit( &k0, &k1, ss.filter_params );

    gsl_integration_qag( &F, k0, k1,
            GSL_INTE_ERR_ABS, GSL_INTE_ERR_REL, GSL_INTE_WS_LEN,
            GSL_INTE_KEY, inte_ws, &r, &err );

    return SQR(ss.norm) * r;
    */

    return 0;

}

void init_ps() {

    double R, sigma8;

    writelog( "init ps norm\n" );

    R = 8.0 / All.HubbleParam;

    global_ss.P = &PowerSpec_Efstathiou;
    global_ss.filter = &top_hat_filter;
    global_ss.filter_k_limit = &top_hat_filter_k_limit;
    global_ss.filter_params = &R;

    /*
    sigma8 = sigma( global_ss );
    All.SigmaNorm = All.Sigma8 / sigma8;
    global_ss.norm = All.SigmaNorm;

    global_ss.dsigma2dmdk = &top_hat_dsigma2dmdk;
    global_ss.dsigma2dmdk_params = &global_ss;

    writelog( "Sigma8: %g, SigmaNorm: %g\n", All.Sigma8, All.SigmaNorm );
    */

}

double dNdM( double a, double M ) {

    /*
    double sigma, dsig2dm, growth, Deltac;

    growth = growth_factor0( a );
    sigma = sigmaM( M ) * growth;
    dsig2dm = dsigma2dm( M ) * ( SQR(growth) / (2*sigma) );

    Deltac = 1.68;

    return ( -All.RhoM/M ) * sqrt( 2.0/M_PI ) * ( Deltac / ( SQR(sigma) ) ) *
        dsig2dm * exp( -SQR(Deltac)/(2*SQR(sigma)) );
        */

    return 0;

}

#define writelog1( a, b ) writelog( "%-35s: %g\n", a, b )
void compute_cosmo_quantities() {
    writelog( "compute cosmology quantities ...\n" );

    All.Time = header.time;
    All.Time2 = SQR( All.Time );
    All.Time3 = CUBE( All.Time );
    All.Hubble_a = hubble_function( All.Time );
    All.RhoCrit = 3 * SQR( All.Hubble_a ) / ( 8*PI*All.G );
    All.RhoBaryon = All.OmegaBaryon * All.RhoCrit;
    All.RhoM = All.Omega0 * All.RhoCrit;

    if ( All.Time > 1e-10 ){
        All.ComDis = comoving_distance( All.Time );
        All.AngDis = angular_distance( All.Time );
        All.LumDis = luminosity_distance( All.Time );
    }
    else
        All.ComDis = All.AngDis = All.LumDis = 0;

    writelog1( "Time", All.Time );
    writelog1( "Time2", All.Time2 );
    writelog1( "Time3", All.Time3 );
    writelog1( "Comoving Distance", All.ComDis );
    writelog1( "Angular Distance", All.AngDis );
    writelog1( "Luminosity Distance", All.LumDis );
    writelog1( "Hubble_a", All.Hubble_a );
    writelog1( "RhoBaryon", All.RhoBaryon );
    writelog1( "RhoCrit", All.RhoCrit );
    writelog1( "RhoM", All.RhoM );

    writelog( "compute cosmology quantities ... done.\n" );
    put_block_line;
}

