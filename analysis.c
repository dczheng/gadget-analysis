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
#define c_in_cgs 2.9979e10
#define HBAR ( 1.05457e-27 * cm * cm * g / s )
#define e2  ( HBAR * LightSpeed / 137.04 )
#define statcoul sqrt( erg * cm )
#define mpc2 ( m_p * LightSpeed * LightSpeed )
#define mec2 ( m_e * LightSpeed * LightSpeed )
#define c2   ( LightSpeed * LightSpeed )

#define GSL_INTE_WS_LEN 1000
#define GSL_INTE_ERR_ABS 0.0
#define GSL_INTE_ERR_REL 1e-3
#define GSL_INTE_KEY GSL_INTEG_GAUSS15

#define SQR(X) ( X*X )
#define CUBE(X) ( X*X*X )

int cpn;
double *cp, *red, *green, *blue, affine[6];
char cb_s, cb_label[100], title[100], xlabel[100], ylabel[100];
static gsl_integration_workspace *inte_ws;

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

void hg_electrons_analysis() {
    double *rho_n, x, y, dx, dy, rho_n_max, rho_n_min,
           log_rho_n_min, log_rho_n_max;
    char buf[100];
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

    if ( access( "./hge_n/", 0 ) == -1 ){
        printf( "create directory ./hge_n/.\n" );
        if ( mkdir( "./hge_n/", 0755) == -1 ){
            printf( "failed create directory ./hge_n/.\n" );
            endrun( 20171130 );
        }
    }
    sprintf( cb_label, "(10^x)" );
    sprintf( buf, "./hge_n/hge_n_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
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
    char buf[200], buf1[100];
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
            "log_rho_max: %g, log_rho_min: %g\n",
            rho_max, rho_min, log_rho_max, log_rho_min );
    for ( i=0; i<PicSize*PicSize; i++ )
        if ( rho[i] > 0 )
            rho[i] = log10( rho[i] );
        else
            rho[i] = log_rho_min - 10 ;

    sprintf( buf1, "./rho_%d", pt );
    if ( access( buf1, 0 ) == -1 ){
        printf( "create directory %s.\n", buf1 );
        if ( mkdir( buf1, 0755) == -1 ){
            printf( "failed create directory %s.\n", buf1 );
            endrun( 20171130 );
        }
    }
    sprintf( cb_label, "(10^x)   g/cm^3" );
    sprintf( buf, "%s/rho_%d_%.2f.png", buf1, pt, RedShift );
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
    char buf[100];
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
        if (SphP[i].MachNumber > mn_max)
            mn_max = SphP[i].MachNumber;
        if (SphP[i].MachNumber < mn_min)
            mn_min = SphP[i].MachNumber;
    }
    printf( "SphP Max MachNumber: %g\n", mn_max );
    printf( "SphP Min MachNumber: %g\n", mn_min );
    mn_max = -DBL_MAX;
    mn_min = DBL_MAX;
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

    if ( access( "./mn/", 0 ) == -1 ){
        printf( "create directory ./mn.\n" );
        if ( mkdir( "./mn/", 0755) == -1 ){
            printf( "failed create directory ./mn.\n" );
            endrun( 20171130 );
        }
    }
    sprintf( cb_label, "(10^x)" );
    sprintf( buf, "./mn/mn_%.2f.png", RedShift );
    giza_open_device( "/png", buf );
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
            "log_rho_max: %g, log_rho_min: %g\n",
            rho_max, rho_min, log_rho_max, log_rho_min );
    for ( i=0; i<PicSize*PicSize; i++ )
        if ( rho[i] > 0 )
            rho[i] = log10( rho[i] );
        else
            rho[i] = log_rho_min - 10 ;

    if ( access( "./gas_rho/", 0 ) == -1 ){
        printf( "create directory ./gas_rho.\n" );
        if ( mkdir( "./gas_rho", 0755) == -1 ){
            printf( "failed create directory ./gas_rho.\n" );
            endrun( 20171130 );
        }
    }
    sprintf( cb_label, "(10^x)   g/cm^3" );
    sprintf( buf, "./gas_rho/gas_rho_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
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

void magnetic_field_analysis() {
    double *mag, mag_max, mag_min, log_mag_max, log_mag_min,
           dx, dy, x, y;
    int i, j, xi, yi;
    char buf[100];
    puts( "magnetic field analysis ...");
    mag = malloc( sizeof( double ) * PicSize * PicSize );
    memset( mag, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = BoxSize / PicSize;
    mag_max = -DBL_MAX;
    mag_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ){
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        mag[ xi*PicSize + yi ] += sqrt( pow( SphP[i].B[0], 2.0 ) +
                pow( SphP[i].B[1],2.0 ) + pow( SphP[i].B[2], 2.0 ) );
    }

    for ( i=0; i<PicSize*PicSize; i++ ) {
        if ( mag[i] > 0 ) {
            if ( mag[i] > mag_max )
                mag_max = mag[i];
            if ( mag[i] < mag_min )
                mag_min = mag[i];
        }
    }
    log_mag_max = log10( mag_max );
    log_mag_min = log10( mag_min );

    printf( "mag_max: %g, mag_min: %g\n"
            "log_mag_max: %g, log_mag_min: %g\n",
            mag_max, mag_min,
            log_mag_max, log_mag_min );

    for ( i=0; i<PicSize*PicSize; i++ ){
        if ( mag[i] > 0 )
            mag[i] = log10( mag[i] );
        else
            mag[i] = log_mag_min - 10;
    }
    if ( access( "./mag/", 0 ) == -1 ){
        printf( "create directory ./mag.\n" );
        if ( mkdir( "./mag", 0755) == -1 ){
            printf( "failed create directory ./mag.\n" );
            endrun( 20171130 );
        }
    }
    sprintf( cb_label, "(10^x)  Gauss" );
    sprintf( buf, "./mag/mag_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, PicSize, 0.0, PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( PicSize, PicSize, mag, 0, PicSize, 0, PicSize, log_mag_min, log_mag_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / MpcFlag );
    sprintf( title, "magnetic field (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, log_mag_min, log_mag_max, cb_label );
    giza_close_device();

    puts( "magnetic field analysis ... done.");
}

double cre_beta_inte( double x, void *params ) {
    double *p = params;
    return pow( x, p[0]-1.0 ) * pow( 1.0-x, p[1]-1.0 );
}

double cre_beta( double a, double b, double x ) {
    double r, err;
    gsl_function F;
    double params[2];
    if ( x<1e-20 )
        return 0;
    if ( b>0.0 ){
        return gsl_sf_beta_inc( a, b, x ) * gsl_sf_beta( a, b );
    }
    //printf( "debug for cre_beta: a=%g, b=%g, x=%g\n", a, b, x );
    F.function = &cre_beta_inte;
    F.params = params;
    params[0] = a;
    params[1] = b;
    gsl_integration_qag( &F, 1e-20, x,
            GSL_INTE_ERR_ABS, GSL_INTE_ERR_REL, GSL_INTE_WS_LEN, GSL_INTE_KEY,
            inte_ws, &r, &err );
    return r;
}

double cre_mean_kinetic_energy( double alpha, double q ) {
    return ( 0.5 * pow( q, alpha-1.0 ) *
            cre_beta( (alpha-2.0)*0.5, (3.0-alpha)*0.5, 1/(1+q*q) ) +
            sqrt( 1+q*q ) - 1.0 ) * mec2;
}

double cre_tau_synchrotron_radiation( double alpha, double q, double B ) {
    /* B: gauss */
    /* per unit electron mass */
    double beta1, ccool, u_b, T;
    //printf( "u_cmb = %g\n", u_cmb );
    /* u_cmb = a_ard * Tcmb0^4 * (z+1)^4 = b^2 / 8 / pi -> b2 ~ (z+1)^4*/
    u_b = SQR(B) / ( 8*M_PI ) * ( erg/CUBE(cm) );
    //printf( "u_b = %g\n", u_b );
    beta1 = cre_beta( ( alpha-2 ) * 0.5, ( 3-alpha ) * 0.5, 1.0 / (1.0+q*q) );
    ccool = 4.0 / (3.0 * m_e*LightSpeed) * (0.665e-24 * cm*cm) * u_b; /* Thomson cross section: 0.665e-24 cm^2 */
    T = cre_mean_kinetic_energy( alpha, q );
    return T / pow( q, alpha-1 ) / mec2 / ( (alpha-1)*ccool*
            ( pow(q,1-alpha)/( alpha-1 ) + pow(q,3-alpha)/(alpha-3) ) );
}

void radio_radiation_analysis() {
    int i;
    double B, dEdt, q, *P;
    puts( "radio analysis ..." );
    for ( i=0; i<N_Gas; i++ ) {
        B = sqrt( pow( SphP[i].B[0], 2.0 ) +
                pow( SphP[i].B[1],2.0 ) + pow( SphP[i].B[2], 2.0 ) );
        q = SphP[i].CRE_Q0 * pow( SphP[i].Density, 0.3333333 );
        //dEdt = cre_tau_synchrotron_radiation( Alpha, q, B );
        SphP[i].P = SphP[i].CRE_E0 / dEdt;
        SphP[i].P /= ( erg/s );
        printf( "%g\n", SphP[i].CRE_C0 );
        //SphP[i].vL = ELECTRONMASS * B / ( 2 * M_PI * ELECTRONCHARGE * c_in_cgs );
        //printf( "%g\n", SphP[i].vL );
    }
    puts( "radio analysis ... done.");
}


void init_analysis() {
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

    inte_ws = gsl_integration_workspace_alloc( GSL_INTE_WS_LEN );
    puts( sep_str );
}

void free_analysis() {
    free( cp );
    free( red );
    free( green );
    free( blue );
    gsl_integration_workspace_free( inte_ws );
}



void gas_analysis(){
    printf( "\n" );
    mach_analysis();
    printf( "\n" );
    hg_electrons_analysis();
    printf( "\n" );
    gas_density_analysis();
    printf( "\n" );
    pos_analysis( 0 );
    printf( "\n" );
    magnetic_field_analysis();
    printf( "\n" );
    magnetic_field_analysis();
    printf( "\n" );
    radio_radiation_analysis();
    fputs( sep_str, stdout );
}

void dm_analysis(){
    printf( "\n" );
    pos_analysis( 1 );
    fputs( sep_str, stdout );
}
