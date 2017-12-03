#include "allvars.h"

#define PROTONMASS   1.6726e-24
#define ELECTRONMASS 9.10953e-28
#define ELECTRONCHARGE  4.8032e-10
#define BOLTZMANN      1.38066e-16
#define cm ( header.HubbleParam / para.UnitLength_in_cm )
#define g  ( header.HubbleParam / para.UnitMass_in_g )
#define s  ( header.HubbleParam / para.UnitTime_in_s )
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
    fprintf( LogFilefd, "high energy electrons analysis...\n" );
    rho_n = ( double* ) malloc( sizeof(double) * para.PicSize * para.PicSize );
    memset( rho_n, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = BoxSize / para.PicSize;
    dy = BoxSize / para.PicSize;

    for ( i=0; i<N_Gas; i++ ) {
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        rho_n[xi*para.PicSize+yi] += SphP[i].CRE_n0 * SphP[i].Density * m_e / ( g/(cm*cm*cm) );
        //printf( "xi: %d, yi: %d, %g\n", xi, yi, rho_n[xi*para.PicSize+yi] );
    }
    rho_n_max = -DBL_MAX;
    rho_n_min = DBL_MAX;
    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( rho_n[i] > 0 ){
            if ( rho_n[i] < rho_n_min )
                rho_n_min = rho_n[i];
            if ( rho_n[i] > rho_n_max )
                rho_n_max = rho_n[i];
        }
    }
    log_rho_n_max = log10( rho_n_max );
    log_rho_n_min = log10( rho_n_min );
    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( rho_n[i] > 0 )
            rho_n[i] = log10( rho_n[i] );
        else
            rho_n[i] = log_rho_n_min - 10;
    }
    fprintf( LogFilefd, "rho_n_max: %g, rho_n_min: %g \n"
             "log_rho_n_max: %g log_rho_n_min: %g\n",
             rho_n_max, rho_n_min, log_rho_n_max, log_rho_n_min );

    if ( ThisTask == 0 )
    if ( access( "./hge_n/", 0 ) == -1 ){
        printf( "create directory ./hge_n/ by task %d\n", ThisTask );
        if ( mkdir( "./hge_n/", 0755) == -1 ){
            printf( "failed create directory ./hge_n/.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    sprintf( cb_label, "(10^x)" );
    sprintf( buf, "./hge_n/hge_n_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, rho_n, 0, para.PicSize, 0, para.PicSize,
            log_rho_n_min, log_rho_n_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "high energy electron number densit (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, log_rho_n_min, log_rho_n_max, cb_label );
    giza_close_device();

    free( rho_n );
    fprintf( LogFilefd, "high energy electrons analysis...done.\n" );
}

void pos_analysis( int pt ){
    FILE *fd;
    int i,j, xi, yi;
    double *rho, x, y, dx, dy, rho_max, rho_min,
           log_rho_max, log_rho_min;
    char buf[200], buf1[100];
    long num, offset;
    fprintf( LogFilefd, "particle %d positin analysis ...\n", pt );
    rho = ( double* ) malloc( sizeof(double) * para.PicSize * para.PicSize );
    memset( rho, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
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
            rho[ xi*para.PicSize + yi ] += P[offset+i].Mass / (dx*dy*BoxSize) / ( g/(cm*cm*cm) );
        else
            rho[ xi*para.PicSize + yi ] += header.mass[pt] / (dx*dy*BoxSize) / ( g/(cm*cm*cm) );
    }
    for ( i=0; i<para.PicSize*para.PicSize; i++ ) {
        if ( rho[i] > 0 ) {
            if ( rho[i] > rho_max )
                rho_max = rho[i];
            if ( rho[i] < rho_min )
                rho_min = rho[i];
        }
    }
    log_rho_max = log10( rho_max );
    log_rho_min = log10( rho_min );

    fprintf( LogFilefd, "rho_max: %g, rho_min: %g\n"
            "log_rho_max: %g, log_rho_min: %g\n",
            rho_max, rho_min, log_rho_max, log_rho_min );
    for ( i=0; i<para.PicSize*para.PicSize; i++ )
        if ( rho[i] > 0 )
            rho[i] = log10( rho[i] );
        else
            rho[i] = log_rho_min - 10 ;

    sprintf( buf1, "./rho_%d", pt );
    if ( ThisTask == 0 )
    if ( access( buf1, 0 ) == -1 ){
        printf( "create directory %s by task %d\n", buf1, ThisTask );
        if ( mkdir( buf1, 0755) == -1 ){
            printf( "failed create directory %s.\n", buf1 );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    sprintf( cb_label, "(10^x)   g/cm^3" );
    sprintf( buf, "%s/rho_%d_%.2f.png", buf1, pt, RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, rho, 0, para.PicSize, 0, para.PicSize, log_rho_min, log_rho_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "particle %d density (z=%.2f)", pt, RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, log_rho_min, log_rho_max, cb_label );
    giza_close_device();


    free( rho );
    fprintf( LogFilefd, "particle %d positin analysis ... done.\n", pt );
}

void mach_analysis(){
    int i,j, xi, yi, pt;
    double *mn, x, y, dx, dy, mn_max, mn_min, log_mn_min, log_mn_max;
    char buf[100];
    pt = 0;
    fprintf( LogFilefd, "mach number analysis ...\n" );
    mn = ( double* ) malloc( sizeof(double) * para.PicSize * para.PicSize );
    memset( mn, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
    mn_max = -DBL_MAX;
    mn_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        mn[ xi*para.PicSize + yi ] += SphP[i].MachNumber;
        if (SphP[i].MachNumber > mn_max)
            mn_max = SphP[i].MachNumber;
        if (SphP[i].MachNumber < mn_min)
            mn_min = SphP[i].MachNumber;
    }
    fprintf( LogFilefd, "SphP Max MachNumber: %g\n", mn_max );
    fprintf( LogFilefd, "SphP Min MachNumber: %g\n", mn_min );
    mn_max = -DBL_MAX;
    mn_min = DBL_MAX;
    for ( i=0; i<para.PicSize*para.PicSize; i++ ) {
        if ( mn[i] > 0 ){
            if ( mn[i] > mn_max )
                mn_max = mn[i];
            if ( mn[i] < mn_min )
                mn_min = mn[i];
        }
    }
    log_mn_min = log10( mn_min );
    log_mn_max = log10( mn_max );
    fprintf( LogFilefd, "mn_max: %g, mn_min: %g\n"
            "log_mn_max: %g, log_mn_min: %g\n",
            mn_max, mn_min, log_mn_max, log_mn_min );
    for ( i=0; i<para.PicSize*para.PicSize; i++ )
        if ( mn[i] > 0 )
            mn[i] = log10( mn[i] );
        else
            mn[i] = log_mn_min - 10 ;

    if ( ThisTask == 0 )
    if ( access( "./mn/", 0 ) == -1 ){
        printf( "create directory ./mn by task %d\n", ThisTask );
        if ( mkdir( "./mn/", 0755) == -1 ){
            printf( "failed create directory ./mn.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    sprintf( cb_label, "(10^x)" );
    sprintf( buf, "./mn/mn_%.2f.png", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, mn, 0, para.PicSize, 0, para.PicSize, log_mn_min, log_mn_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "mach number (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, log_mn_min, log_mn_max, cb_label );
    giza_close_device();

    free( mn );
    fprintf( LogFilefd, "mach number analysis ...done.\n" );
}

void gas_density_analysis(){
    int i,j, xi, yi;
    double *rho, x, y, dx, dy, rho_max, rho_min;
    double log_rho_max, log_rho_min;
    char buf[200];
    fprintf( LogFilefd, "gas density analysis ...\n");
    rho = ( double* ) malloc( sizeof(double) * para.PicSize * para.PicSize );
    memset( rho, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
    rho_max = -DBL_MAX;
    rho_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        rho[ xi*para.PicSize + yi ] += SphP[i].Density / ( g/(cm*cm*cm) );
    }
    for ( i=0; i<para.PicSize*para.PicSize; i++ ) {
        if ( rho[i] > 0 ) {
            if ( rho[i] > rho_max )
                rho_max = rho[i];
            if ( rho[i] < rho_min )
                rho_min = rho[i];
        }
    }
    log_rho_max = log10( rho_max );
    log_rho_min = log10( rho_min );

    fprintf( LogFilefd, "rho_max: %g, rho_min: %g\n"
            "log_rho_max: %g, log_rho_min: %g\n",
            rho_max, rho_min, log_rho_max, log_rho_min );
    for ( i=0; i<para.PicSize*para.PicSize; i++ )
        if ( rho[i] > 0 )
            rho[i] = log10( rho[i] );
        else
            rho[i] = log_rho_min - 10 ;

    if ( ThisTask == 0 )
    if ( access( "./gas_rho/", 0 ) == -1 ){
        printf( "create directory ./gas_rho by task %d\n", ThisTask );
        if ( mkdir( "./gas_rho", 0755) == -1 ){
            printf( "failed create directory ./gas_rho.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    sprintf( cb_label, "(10^x)   g/cm^3" );
    sprintf( buf, "./gas_rho/gas_rho_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, rho, 0, para.PicSize, 0, para.PicSize, log_rho_min, log_rho_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "gas density (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, log_rho_min, log_rho_max, cb_label );
    giza_close_device();

    free( rho );
    fprintf( LogFilefd, "gas density analysis ...done.\n");
}

void magnetic_field_analysis() {
    double *mag, mag_max, mag_min, log_mag_max, log_mag_min,
           dx, dy, x, y;
    int i, j, xi, yi;
    char buf[100];
    fprintf( LogFilefd, "magnetic field analysis ...\n" );
    mag = malloc( sizeof( double ) * para.PicSize * para.PicSize );
    memset( mag, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
    mag_max = -DBL_MAX;
    mag_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ){
        x = P[i].Pos[0];
        y = P[i].Pos[1];
        xi = x / dx;
        yi = y / dy;
        mag[ xi*para.PicSize + yi ] += sqrt( pow( SphP[i].B[0], 2.0 ) +
                pow( SphP[i].B[1],2.0 ) + pow( SphP[i].B[2], 2.0 ) );
    }

    for ( i=0; i<para.PicSize*para.PicSize; i++ ) {
        if ( mag[i] > 0 ) {
            if ( mag[i] > mag_max )
                mag_max = mag[i];
            if ( mag[i] < mag_min )
                mag_min = mag[i];
        }
    }
    log_mag_max = log10( mag_max );
    log_mag_min = log10( mag_min );

    fprintf( LogFilefd, "mag_max: %g, mag_min: %g\n"
            "log_mag_max: %g, log_mag_min: %g\n",
            mag_max, mag_min,
            log_mag_max, log_mag_min );

    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( mag[i] > 0 )
            mag[i] = log10( mag[i] );
        else
            mag[i] = log_mag_min - 10;
    }
    if ( ThisTask == 0 )
    if ( access( "./mag/", 0 ) == -1 ){
        printf( "create directory ./mag by task %d\n", ThisTask);
        if ( mkdir( "./mag", 0755) == -1 ){
            printf( "failed create directory ./mag.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    sprintf( cb_label, "(10^x)  Gauss" );
    sprintf( buf, "./mag/mag_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, mag, 0, para.PicSize, 0, para.PicSize, log_mag_min, log_mag_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "magnetic field (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, log_mag_min, log_mag_max, cb_label );
    giza_close_device();

    fprintf( LogFilefd, "magnetic field analysis ... done.\n" );
}

double cre_beta_inte( double x, void *params ) {
    double *p = params;
    //printf( "%g %g\n", p[0], p[1] );
    return pow( x, p[0]-1.0 ) * pow( 1.0-x, p[1]-1.0 );
}

double cre_beta( double a, double b, double x ) {
    double r, err;
    gsl_function F;
    double params[2];
    //printf( "debug for cre_beta: a=%g, b=%g, x=%g\n", a, b, x );
    if ( x<1e-20 )
        return 0;
    if ( b>0.0 ){
        return gsl_sf_beta_inc( a, b, x ) * gsl_sf_beta( a, b );
    }
    F.function = &cre_beta_inte;
    F.params = params;
    params[0] = a;
    params[1] = b;
    if ( 1.0 - x < 1e-10 )
        x = 1 - 1e-10;
    gsl_integration_qag( &F, 1e-10, x,
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
    int i, j, xi, yi;
    double B, dEdt, q, *p, dx, dy, p_max, p_min,
           log_p_max, log_p_min, x, y;
    char buf[100];
    fprintf( LogFilefd, "radio analysis ...\n" );
    fprintf( LogFilefd, "Alpha =%g\n", para.Alpha );
    p = malloc( sizeof( double ) * para.PicSize * para.PicSize );
    memset( p, 0, sizeof( double ) * para.PicSize * para.PicSize );
    dx = dy = BoxSize / para.PicSize;
    p_max = -DBL_MAX;
    p_min =  DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        if ( SphP[i].CRE_C0 != 0 ) {
            //printf( "%g %g\n", SphP[i].CRE_C0, SphP[i].CRE_Q0 );
            B = sqrt( pow( SphP[i].B[0], 2.0 ) +
                pow( SphP[i].B[1],2.0 ) + pow( SphP[i].B[2], 2.0 ) );
            q = SphP[i].CRE_Q0 * pow( SphP[i].Density, 0.3333333 );
            //printf( "q=%g, B=%g\n", q, B );
            dEdt = cre_tau_synchrotron_radiation( para.Alpha, q, B );
            SphP[i].P = SphP[i].CRE_E0 / dEdt;
            SphP[i].P /= ( erg/s );
            x = P[i].Pos[0];
            y = P[i].Pos[1];
            xi = x / dx;
            yi = y / dy;
            p[ xi*para.PicSize + yi ] += sqrt( pow( SphP[i].B[0], 2.0 ) +
                pow( SphP[i].B[1],2.0 ) + pow( SphP[i].B[2], 2.0 ) );
            //printf( "%g\n", SphP[i].CRE_C0 );
        //SphP[i].vL = ELECTRONMASS * B / ( 2 * M_PI * ELECTRONCHARGE * c_in_cgs );
        //printf( "%g\n", SphP[i].vL );
        }
    }
    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( p[i] > 0 ) {
            if ( p[i] > p_max )
                p_max = p[i];
            if ( p[i] < p_min )
                p_min = p[i];
        }
    }
    log_p_max = log10( p_max );
    log_p_min = log10( p_min );
    fprintf( LogFilefd, "p_max: %g, p_min: %g\n",
            "log_p_max: %g, log_p_min: %g\n",
            p_max, p_min,
            log_p_max, log_p_min );
    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( p[i] > 0 )
            p[i] = log10( p[i] );
        else
            p[i] = log_p_min - 10;
    }
    if ( ThisTask == 0 )
    if ( access( "./rad/", 0 ) == -1 ){
        printf( "create directory ./rad by task %d\n", ThisTask );
        if ( mkdir( "./rad", 0755) == -1 ){
            printf( "failed create directory ./rad.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    sprintf( cb_label, "(10^x)" );
    sprintf( buf, "./rad/rad_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, p, 0, para.PicSize, 0, para.PicSize,
            log_p_min, log_p_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "radio radiation (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, log_p_min, log_p_max, cb_label );
    giza_close_device();
    fprintf( LogFilefd, "radio analysis ... done.\n" );
    free( p );
}

void init_analysis() {
    int i;
    fprintf( LogFilefd, "initialize plot...\n" );
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
    fprintf( LogFilefd, sep_str );
}

void free_analysis() {
    free( cp );
    free( red );
    free( green );
    free( blue );
    gsl_integration_workspace_free( inte_ws );
}

void gas_analysis(){
    mach_analysis();
    fprintf( LogFilefd, "\n" );
    hg_electrons_analysis();
    fprintf( LogFilefd, "\n" );
    gas_density_analysis();
    fprintf( LogFilefd, "\n" );
    pos_analysis( 0 );
    fprintf( LogFilefd, "\n" );
    magnetic_field_analysis();
    fprintf( LogFilefd, "\n" );
    radio_radiation_analysis();
    fprintf( LogFilefd, sep_str );;
}

void dm_analysis(){
    fprintf( LogFilefd, "\n" );
    pos_analysis( 1 );
    fputs( sep_str, LogFilefd );
}
