#include "allvars.h"

void endrun( int ierr ) {
    fprintf( stderr, "EXIT CODE: %i\n", ierr );
    MPI_Abort( MPI_COMM_WORLD, ierr );
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
    if ( MpcFlag != 1 ) {
        MpcFlag = 1000;
    }
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
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &ThisTask );
    MPI_Comm_size( MPI_COMM_WORLD, &NumTask );
    if ( ThisTask == 0 )
        if ( argc < 2 ) {
            fprintf( stderr, "Parameter file is required on command line!\n " );
            endrun( 1 );
        }
#ifdef DEBUG
    signal( SIGSEGV, signal_hander );
#endif
    if ( ThisTask == 0 ){
        time1 = time( NULL );
        tb = localtime( &time1 );
        fprintf( stdout, "Start At: %s", asctime(tb) );
    }

    if ( access( "./gadget-analysis.log/", 0 ) == -1 ) {
        printf( "create directory ./gadget-analysis.log/ by task %d\n", ThisTask );
        if ( mkdir( "./gadget-analysis.log/", 0755 ) == -1 ) {
            printf( "failed create directory ./gadget-analysis.log/\n" );
            endrun( 20171203 );
        }
    }
    sprintf( LogFile, "./gadget-analysis.log/gadget-analysis-%03d.log", ThisTask );
    LogFilefd = fopen( LogFile, "r" );
    init_sep_str();
    read_parameters( argv[1] );
    set_units();
    read_snapshot();
    MPI_Barrier( MPI_COMM_WORLD );
    if ( ThisTask == 0 ) {
        printf( "Read data conpleted on all task.\n",
                "Start analysis ...\n" );
    }
    init_analysis();
    //group_analysis();
    /******************analysis***********************/
    gas_analysis();
    dm_analysis();
    /*************************************************/
    MPI_Barrier( MPI_COMM_WORLD );
    if ( ThisTask == 0 ) {
        printf( "analysis completed on all task.\n",
                "Star free memory ...\n" );
    }
    free_analysis();
    free_memory();
    MPI_Barrier( MPI_COMM_WORLD );
    if ( ThisTask == 0 ) {
        printf( "free memory completed. \n" );
        time2 = time( NULL );
        tb = localtime( &time2 );
        fprintf( stdout, "End At: %s", asctime(tb) );
        fprintf( stdout, "Total Time %i\n", time2-time1 );
    }
    fclose( LogFilefd );
    MPI_Finalize();
}

