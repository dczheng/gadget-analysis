#include "allvars.h"

void init_analysis() {
    int i;
    print_log( "initialize plot..." );
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
        red[i] =   pow( i, 0.9 ) / pow( cpn, 0.9 );
        green[i] =  pow( i, 2 ) / pow( cpn, 2 );
        blue[i] =    pow( i, 0.5 ) / pow( cpn, 0.5 );
    }

    inte_ws = gsl_integration_workspace_alloc( GSL_INTE_WS_LEN );
    print_log( "initialize plot... done." );
    print_log( sep_str );
}

void free_analysis() {
    print_log( "free plot ..." );
    free( cp );
    free( red );
    free( green );
    free( blue );
    gsl_integration_workspace_free( inte_ws );
    print_log( "free plot ... done." );
    print_log( sep_str );
}

void analysis(){
    slice();
    init_projection();
    init_analysis();
    init_kernel_matrix();
    magnetic_field_analysis();
    mach_analysis();
    //fprintf( LogFilefd, "\n" );
    hg_electrons_analysis();
    //fprintf( LogFilefd, "\n" );
    gas_density_analysis();
    pos_analysis( 0 );
    //fprintf( LogFilefd, "\n" );
    //fprintf( LogFilefd, "\n" );
    //fprintf( LogFilefd, "\n" );
    //fprintf( LogFilefd, "\n" );
    //crc_analysis();
    //fprintf( LogFilefd, "\n" );
    //radio_radiation_analysis();
    //fprintf( LogFilefd, sep_str );;
    /*
    fprintf( LogFilefd, "\n" );
    */
    pos_analysis( 1 );
    free_analysis();
    free_kernel_matrix();
}
