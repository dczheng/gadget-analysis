#include "allvars.h"

#define PROTONMASS   1.6726e-24
#define ELECTRONMASS 9.10953e-28
#define ELECTRONCHARGE  4.8032e-10
#define BOLTZMANN      1.38066e-16
#define cm ( header.HubbleParam / UnitLength_in_cm )
#define g  ( header.HubbleParam / UnitMass_in_g )
#define s  ( header.HubbleParam / UnitTime_in_s )
#define erg (g*cm*cm/(s*s))
#define keV (1.602e-9*erg)
#define deg 1.0
#define m_p (PROTONMASS * g)
#define m_e (ELECTRONMASS * g)
#define k_B (BOLTZMANN * erg / deg)
#define LightSpeed (2.9979e10*cm/s)
#define HBAR ( 1.05457e-27 * cm * cm * g / s )
#define e2  ( HBAR * LightSpeed / 137.04 )
#define statcoul sqrt( erg * cm )
#define mpc2 ( m_p * LightSpeed * LightSpeed )
#define mec2 ( m_e * LightSpeed * LightSpeed )
#define c2   ( LightSpeed * LightSpeed )

void magnetic_field_analysis() {
     long i, index;
     double sum, b, bmax, min, max;
     sum = 0;
     index = 0;
     bmax = 0;
     for ( i=0; i<Particle[0].num; i++ ) {
         /*
         if( Particle[0].mag[i*3+0] != 0 ||
             Particle[0].mag[i*3+1] != 0 ||
             Particle[0].mag[i*3+2] != 0 )
         fprintf( stdout, "( %.2f, %.2f, %.2f ): %e %e %e\n",
                 Particle[0].pos[i*3+0],
                 Particle[0].pos[i*3+1],
                 Particle[0].pos[i*3+2],
                 Particle[0].mag[i*3+0],
                 Particle[0].mag[i*3+1],
                 Particle[0].mag[i*3+2] );
                 */
         b = pow( Particle[0].mag[i*3+0], 2 ) +
             pow( Particle[0].mag[i*3+1], 2 ) +
             pow( Particle[0].mag[i*3+2], 2 );
         if ( b>bmax ) {
             index =i;
             bmax = b;
         }
         sum += b;
     }
    fprintf( stdout, "max magnetic field: ( %.2f, %.2f, %.2f ): %e %e %e %e\n",
            Particle[0].pos[index*3+0],
            Particle[0].pos[index*3+1],
            Particle[0].pos[index*3+2],
            Particle[0].mag[index*3+0],
            Particle[0].mag[index*3+1],
            Particle[0].mag[index*3+2], sqrt( bmax ) );
     fprintf( stdout, "Total Magnetic Field: %e\n", sqrt( sum / Particle[0].num ) );
}

int compare_for_sort_group_by_mass( const void *a, const void *b ) {
    return (( ( struct group_struct* )a )->Mass < ( ( struct group_struct* )b )->Mass ) ? 1 : -1;
}

void sort_group_by_mass() {
    qsort( (void*)group, TotNgroups, sizeof( struct group_struct ), &compare_for_sort_group_by_mass );
}

void show_group_info() {
    int i;
    for ( i=0; i<TotNgroups; i++ ) {
        if ( group[i].Mass >= 1e3 )
        fprintf( stdout, "mass: %10.3f, cm: ( %10.3f, %10.3f, %10.3f ), vel: ( %10.3f, %10.3f, %10.3f )\n",
                group[i].Mass, group[i].CM[0], group[i].CM[1], group[i].CM[2],
                group[i].Vel[0], group[i].Vel[1], group[i].Vel[2] );
    }
}

void group_analysis() {
    read_group();
    sort_group_by_mass();
    show_group_info();
    free_group();
}

void velocity_analysis() {
    int pt;
    long i;
    float vmax[3], v;
    pt = 0;
    vmax[0] = Particle[pt].vel[0];
    vmax[1] = Particle[pt].vel[1];
    vmax[2] = Particle[pt].vel[2];
    v = sqrt( pow( vmax[0], 2 ) + pow( vmax[1], 2 ) + pow( vmax[2], 2 ) );
    for ( i=0; i<Particle[pt].num; i++ ) {
        if ( sqrt( pow( Particle[pt].vel[ i*3+0 ], 2 ) +
                   pow( Particle[pt].vel[ i*3+1 ], 2 ) +
                   pow( Particle[pt].vel[ i*3+2 ], 2 ) ) > v ) {
            vmax[0] = Particle[pt].vel[ i*3+0 ];
            vmax[1] = Particle[pt].vel[ i*3+1 ];
            vmax[2] = Particle[pt].vel[ i*3+2 ];
            v = sqrt( pow( vmax[0], 2 ) +
                    pow( vmax[1], 2 ) +
                    pow( vmax[2], 2 ) );
        }
    }
    fprintf( stdout, "Vmax: %f ( %f, %f, %f )\n", v, vmax[0],
            vmax[1], vmax[2] );
}

double gamma_integrand( double t, void *params ) {
    double *x;
    x = ( double* ) params;
    return pow( t, *x-1 ) * exp( -t );
}

double gamma( double x ) {
    double epsabs, epsrel, abserr, result;
    size_t subinter;
    gsl_function F;
    epsabs = epsrel = 1e-8;
    subinter = 100000;
    gsl_integration_workspace *integration_workspace =
        gsl_integration_workspace_alloc( subinter );
    F.function = gamma_integrand;
    F.params = &x;
    gsl_integration_qagiu( &F, 0, epsabs, epsrel, subinter,
            integration_workspace, &result, &abserr );
    gsl_integration_workspace_free( integration_workspace );
    return result;
}

void analysis_radio() {
    /*
    double alpha_e, alpha, sigma_pp, sigma_t, me, mp, GeV, C, e;
    double Bcmb, alpha_v, XHe, Ecmb, Ce, Ae, Bc, C_phy, Eb, v, gamma_fac, q_phy;
    long i;
    alpha = 2.5;
    sigma_pp = 32 * ( 0.96 + exp( 4.4 - 2.4 * alpha ) ) * 1e-28 * 1e4; // cm^2
    me = 9.10938356e-29; // g
    mp = 1.672621898e-24; // g
    sigma_t = 6.6524587158e-29 * 1e4; // cm^2;
    C = 29979245800.0; // m/s
    GeV = 1.602176565e-19 * 1e9 * 1e7; // erg
    Bcmb = 3.24 * pow( 1+redshift, 2 ) * 1e-6; // G
    alpha_e = alpha + 1;
    alpha_v = alpha / 2.0;
    XHe = 0.24;
    e = 1.6021766208e-19; // C
    Ecmb = Bcmb * Bcmb / 8.0 / M_PI;
    v = 150000;
    gamma_fac = gamma( (3*alpha_e-1) / 12.0 ) *
            gamma( (3*alpha_e+7) / 12.0 ) *
            gamma( (alpha_e+5) / 4.0 ) /
            gamma( (alpha_e+7) / 4.0 );
    fprintf( stdout, "alpha = %f\nalpha_e = %f\nalpha_v = %e\n"
            "sigma_pp = %e\nsigma_t = %e\nC = %e\nGeV = %e\n"
            "me = %e\nmp = %e\nBcmb = %e\nXHe = %f\ne = %e\n"
            "Ecmb = %e\ngamma_fac = %e\n",
            alpha, alpha_e, alpha_v, sigma_pp, sigma_t, C, GeV, me, mp,
            Bcmb, XHe, e, Ecmb, gamma_fac );
    Bc = 2.0 * M_PI * pow( me, 3 ) * pow( C, 5 ) * v / ( 3 * e * pow( GeV, 2 ) );
    for ( i=0; i<Particle[0].num; i++ ) {
        if ( Particle[0].c0[i] == 0 ) {
            Particle[0].j[i] = 0;
            continue;
        }
        C_phy =  Particle[0].c0[i] * pow( Particle[0].rho[i], ( alpha - 1 ) * 0.33333 ) *
            Particle[0].rho[i] * 1.989e43 / pow( 3.085678e21, 3) / mp;
        //q_phy = Particle[0].q0[i] * pow( Particle[0].rho[i], 0.3333333 );
        Eb = ( pow( Particle[0].mag[i*3 + 0], 2 ) +
             pow( Particle[0].mag[ i*3 +1 ], 2 ) +
             pow( Particle[0].mag[ i*3 + 2 ], 2 ) ) / 8.0 / M_PI;
        Ce = pow( 16, 2-alpha_e ) / ( alpha_e-2 ) *
            sigma_pp * pow( me,2 ) * pow( C, 4 ) / ( sigma_t*GeV ) *
            C_phy * Particle[0].rho[i] * 1.989e43 / pow( 3.085678e21, 3) *
            Particle[0].elec[i] / ( 1-0.5*0.24 ) / ( Eb + Ecmb ) *
            pow( mp*C*C/GeV, alpha-1 );
        Ae = sqrt( 3.0 * M_PI ) / 32.0 / M_PI * Bc * pow( e, 3 ) / me / pow( C, 2 ) *
            ( alpha_e + 7.0/3.0 ) / ( alpha_e + 1 ) * gamma_fac;
        Particle[0].j[i] = 1e100 * Ae * Ce * pow( Eb/( pow(Bc,2) / 8.0 / M_PI ), ( alpha_v+1 ) / 2.0 );
        fprintf( stdout, "C_phy = %e, q_phy = %e, Eb = %e, Ce = %e, Bc = %e, Ae = %e, j = %e\n",
                C_phy, q_phy, Eb, Ce, Bc, Ae,  Particle[0].j[i] );
    }
    //plot_slice( 0, IO_J );
    */
}

void hg_electrons_analysis() {
    double *rho_n, pos_max[3], pos_min[3], x, y, dx, dy;
    int i,j, pt, xi, yi, zi;
    FILE *fd;
    printf( "high energy electrons analysis...\n" );
    pt = 0;
    rho_n = ( double* ) malloc( sizeof(double) * PicSize * PicSize );
    memset( rho_n, 0, sizeof( double ) * PicSize * PicSize );
    printf( "Particle Number: %ld\n", Particle[pt].num );
    pos_max[0] = pos_max[1] = pos_max[2] = -1e10;
    pos_min[0] = pos_min[1] = pos_min[2] = 1e10;
    dx = header.BoxSize / PicSize;
    dy = header.BoxSize / PicSize;
    for ( i=0; i<Particle[pt].num; i++ ) {
        if ( Particle[pt].pos[i*3+0] > pos_max[0] )
            pos_max[0] = Particle[pt].pos[i*3+0];
        if ( Particle[pt].pos[i*3+1] > pos_max[1] )
            pos_max[1] = Particle[pt].pos[i*3+1];
        if ( Particle[pt].pos[i*3+2] > pos_max[2] )
            pos_max[2] = Particle[pt].pos[i*3+2];
        if ( Particle[pt].pos[i*3+0] < pos_min[0] )
            pos_min[0] = Particle[pt].pos[i*3+0];
        if ( Particle[pt].pos[i*3+1] < pos_min[1] )
            pos_min[1] = Particle[pt].pos[i*3+1];
        if ( Particle[pt].pos[i*3+2] < pos_min[2] )
            pos_min[2] = Particle[pt].pos[i*3+2];
    }
    printf( "pos max: %10g %10g %10g\npos min: %10g %10g %10g\n",
            pos_max[0], pos_max[1], pos_max[2],
            pos_min[0], pos_min[1], pos_min[2] );

    for ( i=0; i<Particle[pt].num; i++ ) {
        x = Particle[pt].pos[i*3+0];
        y = Particle[pt].pos[i*3+1];
        xi = x / dx;
        yi = y / dy;
        rho_n[xi*PicSize+yi] += Particle[pt].cre_n0[i] * Particle[pt].rho[i] * m_e / ( g/(cm*cm*cm) );
        //printf( "xi: %d, yi: %d, %g\n", xi, yi, rho_n[xi*PicSize+yi] );
    }
    fd = fopen( "./hge_n0.txt", "w" );
    for ( i=0; i<PicSize; i++ ){
        for ( j=0; j<PicSize; j++ ) {
            fprintf( fd, "%g ", rho_n[i*PicSize+j] );
        }
        fprintf( fd, "\n" );
    }
    fclose( fd );

    free( rho_n );
    printf( "high energy electrons analysis...done.\n" );
}


void density_analysis( int pt, char *str ){
    FILE *fd;
    int i,j, xi, yi;
    double *rho, x, y, dx, dy, rho_max, rho_min;
    char buf[200];
    printf( "%s density analysis ...\n", str );
    rho = ( double* ) malloc( sizeof(double) * PicSize * PicSize );
    memset( rho, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = header.BoxSize / PicSize;
    rho_max = -1e-10;
    rho_min = 1e10;
    printf( "particle number: %d\n", Particle[pt].num );
    for ( i=0; i<Particle[pt].num; i++ ) {
        x = Particle[pt].pos[i*3+0];
        y = Particle[pt].pos[i*3+1];
        xi = x / dx;
        yi = y / dy;
        rho[ xi*PicSize + yi ] += Particle[pt].rho[i] / ( g/(cm*cm*cm) );
        if ( Particle[pt].rho[i] > rho_max )
            rho_max = Particle[pt].rho[i];
        if ( Particle[pt].rho[i] < rho_min )
            rho_min = Particle[pt].rho[i];
    }
    rho_max /= ( g/(cm*cm*cm) );
    rho_min /= ( g/(cm*cm*cm) );
    printf( "rho_max: %g, rho_min: %g\n", rho_max, rho_min );
    sprintf( buf, "./%s_rho.txt", str );
    fd = fopen( buf, "w" );
    for ( i=0; i<PicSize; i++ ) {
        for ( j=0; j<PicSize; j++ ) {
            fprintf( fd, "%g ", rho[ i*PicSize + j ] );
        }
        fprintf( fd, "\n"  );
    }
    fclose( fd );
    free( rho );
    printf( "%s density analysis ...done.\n", str );
}

void mach_analysis(){
    FILE *fd;
    int i,j, xi, yi, pt;
    double *mn, x, y, dx, dy, mn_max, mn_min;
    pt = 0;
    printf( "mach number analysis ...\n" );
    mn = ( double* ) malloc( sizeof(double) * PicSize * PicSize );
    memset( mn, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = header.BoxSize / PicSize;
    mn_max = -1e-10;
    mn_min = 1e10;
    printf( "particle number: %d\n", Particle[pt].num );
    for ( i=0; i<Particle[pt].num; i++ ) {
        x = Particle[pt].pos[i*3+0];
        y = Particle[pt].pos[i*3+1];
        xi = x / dx;
        yi = y / dy;
        mn[ xi*PicSize + yi ] += Particle[pt].mn[i];
        if ( Particle[pt].mn[i] > mn_max )
            mn_max = Particle[pt].mn[i];
        if ( Particle[pt].mn[i] < mn_min )
            mn_min = Particle[pt].mn[i];
    }
    printf( "mn_max: %g, mn_min: %g\n", mn_max, mn_min );
    fd = fopen( "./mach.txt", "w" );
    for ( i=0; i<PicSize; i++ ) {
        for ( j=0; j<PicSize; j++ ) {
            fprintf( fd, "%g ", mn[ i*PicSize + j ] );
        }
        fprintf( fd, "\n"  );
    }
    fclose( fd );
    free( mn );
    printf( "mach number analysis ...done.\n" );
}


void gas_analysis(){
    printf( "\n" );
    density_analysis( 0, "gas" );
    printf( "\n" );
    hg_electrons_analysis();
    printf( "\n" );
    mach_analysis();
    fputs( sep_str, stdout );
}

void dm_analysis(){
    fputs( sep_str, stdout );
}
