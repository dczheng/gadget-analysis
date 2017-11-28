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

int cpn;
double *cp, *red, *green, *blue, affine[6];
char cb_s, cb_label[100], title[100], xlabel[100], ylabel[100];

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
    double *rho_n, x, y, dx, dy, rho_n_max, rho_n_min,
           log_rho_n_min, log_rho_n_max;
    int i,j, xi, yi, zi;
    printf( "high energy electrons analysis...\n" );
    rho_n = ( double* ) malloc( sizeof(double) * PicSize * PicSize );
    memset( rho_n, 0, sizeof( double ) * PicSize * PicSize );
    dx = BoxSize / PicSize;
    dy = BoxSize / PicSize;

    for ( i=0; i<N_Gas; i++ ) {
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        rho_n[xi*PicSize+yi] += SphP[i].CRE_n0 * SphP[i].Density * m_e / ( g/(cm*cm*cm) );
        //printf( "xi: %d, yi: %d, %g\n", xi, yi, rho_n[xi*PicSize+yi] );
    }
    rho_n_max = -DBL_MAX;
    rho_n_min = DBL_MAX;
    for ( i=0; i<PicSize*PicSize; i++ ){
        if ( rho_n[i] > 0 ){
            if ( rho_n[i] < rho_n_min )
                rho_n_min = rho_n[i];
            if ( rho_n[i] > rho_n_max )
                rho_n_max = rho_n[i];
        }
    }
    log_rho_n_max = log10( rho_n_max );
    log_rho_n_min = log10( rho_n_min );
    for ( i=0; i<PicSize*PicSize; i++ ){
        if ( rho_n[i] > 0 )
            rho_n[i] = log10( rho_n[i] );
        else
            rho_n[i] = log_rho_n_min - 10;
    }
    printf( "rho_n_max: %g, rho_n_min: %g \n"
             "log_rho_n_max: %g log_rho_n_min: %g\n",
             rho_n_max, rho_n_min, log_rho_n_max, log_rho_n_min );

    sprintf( cb_label, "(10^x)" );
    giza_open_device( "/png", "rho_n.png" );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( PicSize, PicSize, rho_n, 0, PicSize, 0, PicSize,
            log_rho_n_min, log_rho_n_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / MpcFlag );
    sprintf( title, "high energy electron number densit (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, log_rho_n_min, log_rho_n_max, cb_label );
    giza_close_device();

    free( rho_n );
    printf( "high energy electrons analysis...done.\n" );
}

void pos_analysis( int pt ){
    FILE *fd;
    int i,j, xi, yi;
    double *rho, x, y, dx, dy, rho_max, rho_min,
           log_rho_max, log_rho_min;
    char buf[200];
    long num, offset;
    printf( "particle %d positin analysis ...\n", pt );
    rho = ( double* ) malloc( sizeof(double) * PicSize * PicSize );
    memset( rho, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = BoxSize / PicSize;
    num = header.npartTotal[pt];
    rho_max = -DBL_MAX;
    rho_min = DBL_MAX;
    num += ( (long long)header.npartTotalHighWord[pt] ) << 32;
    offset = 0;
    for ( i=0; i<pt; i++ ) {
        offset = header.npartTotal[i];
        offset += ( (long long)header.npartTotalHighWord[i] ) << 32;
    }
    for ( i=0; i<num; i++ ) {
        x = P[offset+i].Pos[0];
        y = P[offset+i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        if ( header.mass[pt] == 0 )
            rho[ xi*PicSize + yi ] += P[offset+i].Mass / (dx*dy*BoxSize) / ( g/(cm*cm*cm) );
        else
            rho[ xi*PicSize + yi ] += header.mass[pt] / (dx*dy*BoxSize) / ( g/(cm*cm*cm) );
    }
    for ( i=0; i<PicSize*PicSize; i++ ) {
        if ( rho[i] > 0 ) {
            if ( rho[i] > rho_max )
                rho_max = rho[i];
            if ( rho[i] < rho_min )
                rho_min = rho[i];
        }
    }
    log_rho_max = log10( rho_max );
    log_rho_min = log10( rho_min );

    printf( "rho_max: %g, rho_min: %g\n"
            "log_rho_max: %g, log_rho_min: %g",
            rho_max, rho_min, log_rho_max, log_rho_min );
    for ( i=0; i<PicSize*PicSize; i++ )
        if ( rho[i] > 0 )
            rho[i] = log10( rho[i] );
        else
            rho[i] = log_rho_min - 10 ;

    sprintf( cb_label, "(10^x)   g/cm^3" );
    sprintf( buf, "rho_%d.png", pt );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( PicSize, PicSize, rho, 0, PicSize, 0, PicSize, log_rho_min, log_rho_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / MpcFlag );
    sprintf( title, "particle %d density (z=%.2f)", pt, RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, log_rho_min, log_rho_max, cb_label );
    giza_close_device();


    free( rho );
    printf( "particle %d positin analysis ... done.\n", pt );
}

void mach_analysis(){
    int i,j, xi, yi, pt;
    double *mn, x, y, dx, dy, mn_max, mn_min, log_mn_min, log_mn_max;
    pt = 0;
    printf( "mach number analysis ...\n" );
    mn = ( double* ) malloc( sizeof(double) * PicSize * PicSize );
    memset( mn, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = BoxSize / PicSize;
    mn_max = -DBL_MAX;
    mn_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        mn[ xi*PicSize + yi ] += SphP[i].MachNumber;
    }
    for ( i=0; i<PicSize*PicSize; i++ ) {
        if ( mn[i] > 0 ){
            if ( mn[i] > mn_max )
                mn_max = mn[i];
            if ( mn[i] < mn_min )
                mn_min = mn[i];
        }
    }
    log_mn_min = log10( mn_min );
    log_mn_max = log10( mn_max );
    printf( "mn_max: %g, mn_min: %g\n"
            "log_mn_max: %g, log_mn_min: %g\n",
            mn_max, mn_min, log_mn_max, log_mn_min );
    for ( i=0; i<PicSize*PicSize; i++ )
        if ( mn[i] > 0 )
            mn[i] = log10( mn[i] );
        else
            mn[i] = log_mn_min - 10 ;

    sprintf( cb_label, "(10^x)" );
    giza_open_device( "/png", "mn.png" );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( PicSize, PicSize, mn, 0, PicSize, 0, PicSize, log_mn_min, log_mn_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / MpcFlag );
    sprintf( title, "mach number (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, log_mn_min, log_mn_max, cb_label );
    giza_close_device();

    free( mn );
    printf( "mach number analysis ...done.\n" );
}

void gas_density_analysis(){
    int i,j, xi, yi;
    double *rho, x, y, dx, dy, rho_max, rho_min;
    double log_rho_max, log_rho_min;
    char buf[200];
    printf( "gas density analysis ...\n");
    rho = ( double* ) malloc( sizeof(double) * PicSize * PicSize );
    memset( rho, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = BoxSize / PicSize;
    rho_max = -DBL_MAX;
    rho_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        rho[ xi*PicSize + yi ] += SphP[i].Density / ( g/(cm*cm*cm) );
    }
    for ( i=0; i<PicSize*PicSize; i++ ) {
        if ( rho[i] > 0 ) {
            if ( rho[i] > rho_max )
                rho_max = rho[i];
            if ( rho[i] < rho_min )
                rho_min = rho[i];
        }
    }
    log_rho_max = log10( rho_max );
    log_rho_min = log10( rho_min );

    printf( "rho_max: %g, rho_min: %g\n"
            "log_rho_max: %g, log_rho_min: %g",
            rho_max, rho_min, log_rho_max, log_rho_min );
    for ( i=0; i<PicSize*PicSize; i++ )
        if ( rho[i] > 0 )
            rho[i] = log10( rho[i] );
        else
            rho[i] = log_rho_min - 10 ;

    sprintf( cb_label, "(10^x)   g/cm^3" );
    giza_open_device( "/png", "gas_rho.png" );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( PicSize, PicSize, rho, 0, PicSize, 0, PicSize, log_rho_min, log_rho_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / MpcFlag );
    sprintf( title, "gas density (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, log_rho_min, log_rho_max, cb_label );
    giza_close_device();

    free( rho );
    printf( "gas density analysis ...done.\n");
}

void init_plot() {
    int i;
    puts( "initialize plot..." );
    affine[0] = 1;
    affine[1] = 0;
    affine[2] = 0;
    affine[3] = 1;
    affine[4] = 0;
    affine[5] = 0;
    cb_s = 'R';
    cpn = 30;
    cp = malloc( sizeof(double) * cpn );
    red = malloc( sizeof(double) * cpn );
    green = malloc( sizeof(double) * cpn );
    blue = malloc( sizeof(double) * cpn );
    for ( i=0; i<cpn; i++ ) {
        cp[i] = i / (double)cpn;
        red[i] =   pow( i, 1.5 ) / pow( cpn, 1.5 );
        green[i] =  pow( i, 3 ) / pow( cpn, 33 );
        blue[i] =    pow( i, 1 ) / pow( cpn, 1 );
    }
    puts( sep_str );
}

void free_plot() {
    free( cp );
    free( red );
    free( green );
    free( blue );
}



void gas_analysis(){
    printf( "\n" );
    gas_density_analysis();
    printf( "\n" );
    hg_electrons_analysis();
    printf( "\n" );
    mach_analysis();
    printf( "\n" );
    pos_analysis( 0 );
    fputs( sep_str, stdout );
}

void dm_analysis(){
    printf( "\n" );
    pos_analysis( 1 );
    fputs( sep_str, stdout );
}
