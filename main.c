#include "allvars.h"

void endrun( int ierr ) {
    fprintf( stderr, "EXIT CODE: %i\n", ierr );
    exit( ierr );
}

void set_units() {
    fputs( sep_str, stdout );
    fprintf( stdout, "Set Units... \n" );
    UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    UnitDensity_in_cgs = UnitMass_in_g / pow( UnitLength_in_cm, 3 );
    UnitEnergy_in_cgs = UnitMass_in_g * pow( UnitLength_in_cm,2 ) / pow( UnitTime_in_s, 2 );
    printf( "%-35s: %g\n", "UnitMass_in_g", UnitMass_in_g );
    printf( "%-35s: %g\n", "UnitTime_in_s", UnitTime_in_s );
    printf( "%-35s: %g\n", "UnitLength_in_cm", UnitLength_in_cm );
    printf( "%-35s: %g\n", "UnitDensity_in_cgs", UnitDensity_in_cgs );
    printf( "%-35s: %g\n", "UnitEnergy_in_cgs", UnitEnergy_in_cgs );
    printf( "%-35s: %g\n", "UnitVelocity_in_cm_per_s", UnitVelocity_in_cm_per_s );
    fputs( sep_str, stdout );
}

void init_sep_str() {
    memset( sep_str, '-', SEP_LEN-1 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';
}

int main( int argc, char *argv[] ){
    int i;
    time_t time1, time2;
    struct tm *tb;
    char tmp[100];
    if ( argc < 2 ) {
        fprintf( stderr, "Parameter file is required on command line!\n " );
        endrun( 1 );
    }
#ifdef DEBUG
    signal( SIGSEGV, signal_hander );
#endif
    time1 = time( NULL );
    tb = localtime( &time1 );
    fprintf( stdout, "Start At: %s", asctime(tb) );
    init_sep_str();
    read_parameters( argv[1] );
    set_units();
    read_snapshot();
    //group_analysis();
    /******************analysis***********************/
    gas_analysis();
    dm_analysis();
    /*************************************************/
    time2 = time( NULL );
    tb = localtime( &time2 );
    free_memory();
    fprintf( stdout, "End At: %s", asctime(tb) );
    fprintf( stdout, "Total Time %i\n", time2-time1 );
}

