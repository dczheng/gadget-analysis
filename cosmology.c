#include "allvars.h"

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

double dsigma_dk( double k, void *params ) {

    double p, *R, kR, w;

    R = params;

    p = PowerSpec_Efstathiou( k );

    kR = *R * k;

    if ( kR < 1e-4 )
        w = 1;
    else
        w = 3.0 * ( sin(kR) / CUBE( kR ) - cos( kR ) / SQR(kR) );

    return SQR(k) * p * SQR(w);

}

double sigmaR( double R ) {

    double r, err, k0, k1;

    gsl_function F;
    F.function = dsigma_dk;
    F.params = &R;


    k0 = 1e-99 / R;
    k1 = 350.0 / R;

    gsl_integration_qag( &F, k0, k1,
            GSL_INTE_ERR_ABS, GSL_INTE_ERR_REL, GSL_INTE_WS_LEN,
            GSL_INTE_KEY, inte_ws, &r, &err );

    return All.SigmaNorm * sqrt( r );

}

double sigmaM( double M ) {

    double R;
    R = 3 * M / pow( 3 * M_PI * All.RhoM, 1.0/3.0 );

    return sigmaR( R );

}

void init_sigma_norm() {

    double R, sigma8;

    writelog( "init sigma norm\n" );
    R = 8.0 / All.HubbleParam;

    All.SigmaNorm = 1;

    sigma8 = sigmaR( R );

    All.SigmaNorm = All.Sigma8 / sigma8;

    writelog( "Sigma8: %g, SigmaNorm: %g\n", All.Sigma8, All.SigmaNorm );

}

double dsigma2dmdk( double k, void *params) {

    double p, w, dwdr, kR, *R, kR3, kR2, drdm, R2, *d2fac;

    R = params;
    d2fac = params;
    d2fac ++;

    p = PowerSpec_Efstathiou( k );

    R2 = *R * *R;
    kR = *R * k;
    kR2 = kR * kR;
    kR3 = kR2 * kR;

    if ( kR < 1e-4 )
        w = 1;
    else
        w = 3.0 * ( sin(kR) / kR3 - cos( kR ) / kR2 );

    if ( kR < 1e-10 )
        dwdr = 0;
    else
        dwdr = 9 * cos( kR ) * k / kR3 + 3.0 * sin(kR) * ( 1-3/kR2 ) / (kR* *R);

    drdm = 1 / ( 4 * M_PI * All.RhoM * R2 );

    return SQR(k) * p * 2 * w * dwdr * drdm * *d2fac;

}

double dsigma2dm( double M ) {

    double R, p[2], d2fac, k0, k1, r, err;

    R = 3 * M / pow( 3 * M_PI * All.RhoM, 1.0/3.0 );
    d2fac = M*1e4 / sigmaM( M );

    gsl_function F;
    F.function = &dsigma2dmdk;
    p[0] = R;
    p[1] = d2fac;
    F.params = p;

    k0 = 1e-99 / R;
    k1 = 350.0 / R;

    gsl_integration_qag( &F, k0, k1,
            GSL_INTE_ERR_ABS, GSL_INTE_ERR_REL, GSL_INTE_WS_LEN,
            GSL_INTE_KEY, inte_ws, &r, &err );

    return SQR(All.SigmaNorm) * r / d2fac;


}

double dNdM( double a, double M ) {

    double sigma, dsig2dm, growth, Deltac;

    growth = growth_factor0( a );
    sigma = sigmaM( M ) * growth;
    dsig2dm = dsigma2dm( M ) * ( SQR(growth) / (2*sigma) );

    Deltac = 1.68;

    return ( -All.RhoM/M ) * sqrt( 2.0/M_PI ) * ( Deltac / ( SQR(sigma) ) ) *
        dsig2dm * exp( -SQR(Deltac)/(2*SQR(sigma)) );

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

