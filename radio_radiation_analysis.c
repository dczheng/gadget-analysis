#include "allvars.h"

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
            sqrt( 1+q*q ) - 1.0 ) * cgs_mec2;
}

double cre_tau_synchrotron_radiation( double alpha, double q, double B ) {
    /* B: gauss */
    /* per unit electron mass */
    double beta1, ccool, u_b, T;
    //printf( "u_cmb = %g\n", u_cmb );
    /* u_cmb = a_ard * Tcmb0^4 * (z+1)^4 = b^2 / 8 / pi -> b2 ~ (z+1)^4*/
    u_b = SQR(B) / ( 8*M_PI ) * ( cgs_erg/CUBE(cgs_cm) );
    //printf( "u_b = %g\n", u_b );
    beta1 = cre_beta( ( alpha-2 ) * 0.5, ( 3-alpha ) * 0.5, 1.0 / (1.0+q*q) );
    ccool = 4.0 / (3.0 * cgs_me*cgs_c) * (0.665e-24 * cgs_cm*cgs_cm) * u_b; /* Thomson cross section: 0.665e-24 cm^2 */
    T = cre_mean_kinetic_energy( alpha, q );
    return T / pow( q, alpha-1 ) / cgs_mec2 / ( (alpha-1)*ccool*
            ( pow(q,1-alpha)/( alpha-1 ) + pow(q,3-alpha)/(alpha-3) ) );
}

void radio_radiation_analysis() {
    int i, j, xi, yi;
    double B, dEdt, q, *p, dx, dy, p_max, p_min,
           log_p_max, log_p_min, x, y,
           glob_log_p_max, glob_log_p_min;
    char buf[100];
    print_log( "radio analysis ..." );
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
            SphP[i].P /= ( cgs_erg/cgs_s );
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
    sprintf( LogBuf, "p_max: %g, p_min: %g\n",
            "log_p_max: %g, log_p_min: %g",
            p_max, p_min,
            log_p_max, log_p_min );
    print_log( LogBuf );
    if ( ThisTask == 0 )
    if ( access( "./rad/", 0 ) == -1 ){
        print_log( "create directory `./rad` by task 0" );
        if ( mkdir( "./rad", 0755) == -1 ){
            printf( "failed create directory ./rad.\n" );
            endrun( 20171130 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Reduce( &log_p_max, &glob_log_p_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_p_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Reduce( &log_p_min, &glob_log_p_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Bcast( &glob_log_p_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    sprintf( LogBuf, "glob_log_p_max: %g, glob_log_p_min: %g",
            glob_log_p_max, glob_log_p_min );
    print_log( LogBuf );
    for ( i=0; i<para.PicSize*para.PicSize; i++ ){
        if ( p[i] > 0 )
            p[i] = log10( p[i] );
        else
            p[i] = glob_log_p_min;
    }
    sprintf( cb_label, "(10^x)" );
    sprintf( buf, "./rad/rad_%.2f\n", RedShift );
    giza_open_device( "/png", buf );
    giza_set_environment( 0.0, para.PicSize, 0.0, para.PicSize, 1, -1 );
    giza_set_colour_table( cp, red, green, blue, cpn, 1, 1 );
    giza_render( para.PicSize, para.PicSize, p, 0, para.PicSize, 0, para.PicSize,
            glob_log_p_min, glob_log_p_max, 0, affine );
    sprintf( xlabel, "%g Mpc", BoxSize / para.MpcFlag );
    sprintf( title, "radio radiation (z=%.2f)", RedShift );
    giza_label( xlabel, "", title );
    giza_colour_bar( &cb_s, 1, 3, glob_log_p_min, glob_log_p_max, cb_label );
    fp_tmp = stdout;
    stdout = LogFilefd;
    giza_close_device();
    stdout = fp_tmp;
    free( p );
    print_log( "radio analysis ... done." );
    print_log( sep_str );
}

