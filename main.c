#include "allvars.h"

void endrun( int ierr ) {
    fprintf( stderr, "EXIT CODE: %i\n", ierr );
    MPI_Abort( MPI_COMM_WORLD, ierr );
    exit( ierr );
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

    if ( ThisTask == 0 )
    if ( access( "./gadget-analysis.log/", 0 ) == -1 ) {
        printf( "create directory ./gadget-analysis.log/ by task %d\n", ThisTask );
        if ( mkdir( "./gadget-analysis.log/", 0755 ) == -1 ) {
            printf( "failed create directory ./gadget-analysis.log/\n" );
            endrun( 20171203 );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    sprintf( LogFile, "./gadget-analysis.log/gadget-analysis-%03d.log", ThisTask );
    LogFilefd = fopen( LogFile, "w" );
    init_sep_str();
    read_parameters( argv[1] );
    if ( ThisTask == 0 ){
        printf( sep_str );
        printf( "open log file\n" );
        printf( sep_str );
    }
    set_units();
    read_snapshot();
    MPI_Barrier( MPI_COMM_WORLD );
    if ( ThisTask == 0 ) {
        printf( "Read data completed on all task.\n" );
        printf( sep_str );
        printf( "Start analysis ...\n" );
    }
    init_analysis();
    //group_analysis();
    /******************analysis***********************/
    gas_analysis();
    dm_analysis();
    /*************************************************/
    MPI_Barrier( MPI_COMM_WORLD );
    if ( ThisTask == 0 ) {
        printf( "analysis completed on all task.\n" );
        printf( sep_str );
        printf(  "Star free memory ...\n" );
    }
    free_analysis();
    free_memory();
    MPI_Barrier( MPI_COMM_WORLD );
    if ( ThisTask == 0 ) {
        printf( "free memory completed. \n" );
        time2 = time( NULL );
        tb = localtime( &time2 );
        printf( sep_str );
        fprintf( stdout, "End At: %s", asctime(tb) );
        fprintf( stdout, "Total Time %i\n", time2-time1 );
    }
    fclose( LogFilefd );
    MPI_Finalize();
}

