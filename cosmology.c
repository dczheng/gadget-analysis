#include "allvars.h"

sigma_struct ps_ss;
int init_ps_flag=1;

double E_a ( double a ) {

    double E2;
    E2 =  Omega0 / ( a*a*a ) +
          ( 1-Omega0-OmegaLambda ) / ( a*a ) +
          OmegaLambda;
    return sqrt( E2 );

}

double hubble_function( double a ) {

    return Hubble * E_a( a );

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
    //d *= LIGHT_SPEED / ( HUBBLE * HubbleParam );
    d *= guc.c / ( Hubble );
    //printf( "%g, %g\n", a, d );
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
    return  Omega0 /
        ( a + Omega0 * ( 1-a )  + OmegaLambda * ( CUBE(a)-a ) );

}

double OmegaLambda_a( double a ) {
    // only for Omega_r = Omega_k = 0
    // C.P. Dullemond
    return  OmegaLambda * CUBE(a) /
        ( a + Omega0 * ( 1-a )  + OmegaLambda * ( CUBE(a)-a ) );

}

double growth_factor( double a ){
    // C.P. Dullemond
    //
    double Om, OL;
    Om = OmegaM( a );
    OL = OmegaLambda_a( a );
    return 2.5 * Om / ( pow( Om, 4/7 ) - OL + ( 1+0.5*Om ) * ( 1+1/70.0*OL ) );

}

double growth_factor0( double a ) {

    return growth_factor( a ) / growth_factor( 1 );

}

double PowerSpec_Efstathiou(double k) {

  double AA, BB, CC, nu, ShapeGamma, p;

  //ShapeGamma = 0.25;
  ShapeGamma = 0.21;
  AA = 6.4 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm); // convert unit
  BB = 3.0 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  CC = 1.7 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);

  nu = 1.13;

  p =  k / pow(1 + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);

  return p;

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
    *k1 = 350.0 / *R;


}

double dsigma_dk( double k, void *params ) {

    double p, w;

    sigma_struct *ss;
    ss = params;

    p = (*(ss->P))( k );

    w =(*(ss->filter))( k, ss->FilterParams );

    //printf( "k: %e, p: %e, w: %e\n", k, p, w );

    return SQR(k) * p * SQR(w);

}

double sigma( sigma_struct ss ) {

    double r, err, k0, k1;

    gsl_function F;
    F.function = &dsigma_dk;
    F.params = &ss;
    ss.FilterKLimit( &k0, &k1, ss.FilterParams );

    /*
    for ( k0=1e-10; k0<1e2; k0*=5 )
        printf( "k: %e P: %e\n", k0, PowerSpec_Efstathiou( k0 ) );
    printf( "k0: %g, k1: %g, dsigma_dk: %g, %g\n", k0, k1,
            dsigma_dk( k0, &ss ), dsigma_dk( k1, &ss ) );
            */

    //r = qtrap( &dsigma_dk, &ss, k0, k1, 0.01 );
    gsl_integration_qag( &F, k0, k1,
            GSL_INTE_ERR_ABS, 1e-6, GSL_INTE_WS_LEN,
            GSL_INTE_KEY, inte_ws, &r, &err );

    //printf( "[in sigma] k0: %g, k1: %g, r: %g, err: %g\n", k0, k1, r, err );
    //printf( "dsigma_dk: %g %g\n", dsigma_dk( k0, &ss ), dsigma_dk( k1, &ss ) );

    return sqrt( r ) * ss.norm;

}

void test_sigma() {

    sigma_struct ss;
    double R;

    printf( "test sigma ...\n" );

    R = 8000.0;
    ss.filter = &top_hat_filter;
    ss.FilterParams = &R;
    ss.FilterKLimit = &top_hat_filter_k_limit;
    ss.P = &PowerSpec_Efstathiou;

    printf( "Sigma8: %g\n", sigma( ss ) );

    endrun( 20181025 );

}

double top_hat_dsigma2dmdk( double k, void *params) {

    double p, w, dwdr, kR, *R, kR3, kR2, drdm, R2;

    sigma_struct *ss;
    ss = params;

    R = ss->FilterParams;

    p = (*(ss->P))( k );

    R2 = *R * *R;
    kR = *R * k;
    kR2 = kR * kR;
    kR3 = kR2 * kR;

    w = (*(ss->filter))( k, ss->FilterParams );

    if ( kR < 1e-10 )
        dwdr = 0;
    else
        dwdr = 9 * cos( kR ) * k / kR3 + 3.0 * sin(kR) * ( 1-3/kR2 ) / (kR* *R);

    drdm = 1 / ( 4 * M_PI * RhoM * R2 );

    return SQR(k) * p * 2 * w * dwdr * drdm;

}

double dsigma2dm( sigma_struct ss ) {

    double k0, k1, r, err;

    gsl_function F;
    F.function = ss.dsigma2dmdk;
    F.params = ss.dsigma2dmdk_params;

    (*ss.FilterKLimit)( &k0, &k1, ss.FilterParams );

    gsl_integration_qag( &F, k0, k1,
            GSL_INTE_ERR_ABS, 1e-6, GSL_INTE_WS_LEN,
            GSL_INTE_KEY, inte_ws, &r, &err );

    return SQR(ss.norm) * r;

}

void MtoFilterParams( double M, void *params ) {

    double R, *p;
    sigma_struct *ss;
    ss = params;
    R =  pow( 3.0 * M / (4 * M_PI * RhoM), 1.0/3.0 );
    p = ss->FilterParams;
    p[0] = R;

}

void init_ps() {

    double R, sigma8;

    R = 8000.0;

    writelog( "init ps ...\n" );
    ps_ss.P = &PowerSpec_Efstathiou;
    ps_ss.filter = &top_hat_filter;
    ps_ss.FilterKLimit = &top_hat_filter_k_limit;
    ps_ss.dsigma2dmdk = &top_hat_dsigma2dmdk;
    ps_ss.dsigma2dmdk_params = &ps_ss;
    ps_ss.MtoFilterParams = &MtoFilterParams;

    ps_ss.norm = 1;
    ps_ss.FilterParams = &R;

    sigma8 = sigma( ps_ss );
    ps_ss.norm = All.Sigma8 / sigma8;

    writelog( "Sigma8: %g, sigma: %g, norm: %g\n", All.Sigma8,  sigma8, ps_ss.norm );

}

double PS_dndM( double a, double M ) {

    double p[2];

    if ( init_ps_flag ) {
        init_ps();
        init_ps_flag = 0;
    }

    double sigmaM, dsig2dm, growth, Deltac;

    ps_ss.FilterParams = p;

    (*ps_ss.MtoFilterParams)( M, &ps_ss );

    //printf( "R: %g\n", *(double*)(ps_ss.FilterParams) );

    growth = growth_factor0( a );
    sigmaM = sigma( ps_ss ) * growth;

    dsig2dm =  dsigma2dm( ps_ss ) * ( SQR(growth) / (2*sigmaM) );

    Deltac = 1.68;

    //printf( "growth: %g, sigmaM: %g, dsig2dm: %g\n", growth, sigmaM, dsigma2dm( ps_ss ) * RhoM / M );

    return ( -RhoM/M ) * sqrt( 2.0/M_PI ) * ( Deltac / ( SQR(sigmaM) ) ) *
        dsig2dm * exp( -SQR(Deltac)/(2*SQR(sigmaM)) );

}

void test_ps() {

    printf( "test ps ...\n" );

    double M = 1e3, a = 1, m;


    printf( "RhoM: %g\n", RhoM );
    for( M=1e8; M<1e14; M*=1.1 ) {
        m = M / 1e10;
        printf( "M: %g, dNdM: %g\n", M, PS_dndM( a, m ) * 1e9 / 1e10 * ( CUBE(HubbleParam) * HubbleParam ) );
    }

    endrun( 20181026 );

}

#define writelog1( a, b ) writelog( "%-35s: %g\n", a, b )
void compute_cosmo_quantities() {
    writelog( "compute cosmology quantities ...\n" );

    Time = header.time;
    Time2 = SQR( Time );
    Time3 = CUBE( Time );
    Hubble_a = hubble_function( Time );
    RhoCrit = 3 * SQR( Hubble_a ) / ( 8*PI*G ); 
    RhoBaryon = All.OmegaBaryon / Omega0 * OmegaM(Time) * RhoCrit;
    RhoM = Omega0 * RhoCrit;
    //test_cos();

    if ( Time > 1e-10 ){
        ComDis = comoving_distance( Time );
        AngDis = angular_distance( Time );
        LumDis = luminosity_distance( Time );
    }
    else
        ComDis = AngDis = LumDis = 0;

    writelog1( "Time", Time );
    writelog1( "Time2", Time2 );
    writelog1( "Time3", Time3 );
    writelog1( "Comoving Distance", ComDis );
    writelog1( "Angular Distance", AngDis );
    writelog1( "Luminosity Distance", LumDis );
    writelog1( "Hubble_a", Hubble_a );
    writelog1( "RhoBaryon", RhoBaryon );
    writelog1( "RhoBaryon[cgs]", RhoBaryon * g2c.g / CUBE( g2c.cm ) );
    writelog1( "RhoCrit", RhoCrit );
    writelog1( "RhoCrit[cgs]", RhoCrit * g2c.g / CUBE( g2c.cm ) );
    writelog1( "RhoM", RhoM );

    writelog( "compute cosmology quantities ... done.\n" );
    put_sep;
}

void test_cos() {

    double z, La, Ll, Lc, a;
    HubbleParam = 0.68;
    OmegaLambda = 0.698;
    Omega0 = 0.302;
    set_units();

    if ( ThisTask == 0 ) {
        for ( z=0.1; z<1; z *= 1.1) {
            a = 1 / ( 1+z );
            Lc = comoving_distance( a );
            La = angular_distance( a );
            Ll = luminosity_distance( a );
            printf( "%g, %g, %g, %g, %g, %g, %g\n",
                    z, a, Lc, La, Ll,
                    La / Lc, Ll / Lc );
        }
    }

    do_sync( "" );
    endrun( 20190115 );

}

