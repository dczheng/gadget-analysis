#include "allvars.h"

double particle_df( double p, void *params ) {

    double *pa;
    pa = params;

    if ( p < pa[2] || p > pa[3] )
        return 0;

    return pa[0] * pow( p, -pa[1] );

}

double particle_radio2( double nu,  SphParticleData *part ) {

    double r, params[4], B;

    params[0] = part->CRE_C;
    params[0] = params[0] * part->Density / ( ELECTRON_MASS /(g2c.g) );
    params[0] /= CUBE( g2c.cm );

    params[1] = part->CRE_Alpha;
    params[2] = part->CRE_qmin;
    params[3] = part->CRE_qmax;


    if (  params[ 0] *
            params[1] *
            params[2] *
            params[3] == 0 )
        return 0;

    B = pow( part->B[0], 2 ) + pow( part->B[1], 2 ) + pow( part->B[2], 2 );
    B += SQR( BCMB0 ) * pow( All.Time, -2 );
    B = sqrt( B );

    //printf( "B: %g\n", B );

    r = radio( &particle_df, params, B, nu, params[2], params[3] );

    r = r * ( 4.0/3.0 * PI * CUBE( All.SofteningTable[0] * g2c.cm ) );

    return r;

}

double particle_radio( double nu, long i ) {

    return particle_radio2( nu, &SphP[i] );

}

void test_qmax() {

    SphParticleData part;
    double nu_min, nu_max, nu, dlognu, P, qmax_min, qmax_max, dlogqmax, B;
    int i, N, qmaxn, k;
    FILE *fd;

    part.CRE_C = 1 * CUBE(g2c.cm) * ( ELECTRON_MASS/(g2c.g) );
    part.CRE_Alpha = 2.01;
    part.CRE_qmin = 1;
    part.CRE_qmax = 1e5;
    part.Density = 1;
    B = BCMB0;
    part.B[0] = sqrt( SQR(B) - SQR(BCMB0) * pow( All.Time, -2 ) );
    part.B[1] = part.B[2] = 0;

    N = 100;
    nu_min = 1e6;
    nu_max = 1e9;
    dlognu = log10( nu_max/nu_min ) / ( N-1 );

    qmaxn = 5;
    qmax_min = 5e3;
    qmax_max = 1e5;
    dlogqmax = log10( qmax_max/qmax_min ) / ( qmaxn-1 );

    fd = fopen( "./qmax_test.dat", "w" );

    fprintf( fd, "0 " );

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

        fprintf( fd, "%g ",
                part.CRE_qmax );

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

void test_radio() {


    All.Time = 1;
    All.HubbleParam = 0.7;
    set_units();
    put_block_line;

#ifdef RADIO_F_INTERP
    init_tab_F();
#endif



    if ( ThisTask == 0 ) {
        test_F();
        //test_qmax();
    }

    MPI_Barrier( MPI_COMM_WORLD );
#ifdef RADIO_F_INTERP
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

