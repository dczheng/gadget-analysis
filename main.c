#include "allvars.h"

void endrun( int ierr ) {
    fprintf( stderr, "EXIT CODE: %i\n", ierr );
    MPI_Abort( MPI_COMM_WORLD, ierr );
    exit( ierr );
}

void init_sep_str() {
    memset( sep_str, '-', SEP_LEN-1 );
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
    time1 = time( NULL );
    tb = localtime( &time1 );

    init_sep_str();
    if ( ThisTask == 0 ) {
        printf( "%s\n", sep_str );
        if ( access( "./gadget-analysis.log/", 0 ) == -1 ) {
            printf( "create directory ./gadget-analysis.log/ by task 0\n" );
            if ( mkdir( "./gadget-analysis.log/", 0755 ) == -1 ) {
                printf( "failed create directory ./gadget-analysis.log/\n" );
                endrun( 20171203 );
            }
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    sprintf( LogFile, "./gadget-analysis.log/gadget-analysis-%03d.log", ThisTask );
    LogFilefd = fopen( LogFile, "w" );

    print_log( "open log file" );
    sprintf( LogBuf, "Start At: %s", asctime(tb) );
    LogBuf[strlen( LogBuf )-1] = '\0';
    print_log( LogBuf );

    print_log( sep_str );
    read_parameters( argv[1] );
    read_snapshot();

    set_units();
    //group_analysis();
    /******************analysis***********************/
    analysis();
    MPI_Barrier( MPI_COMM_WORLD );
    /*************************************************/

    free_memory();
    MPI_Barrier( MPI_COMM_WORLD );

    time2 = time( NULL );
    tb = localtime( &time2 );

    sprintf( LogBuf, "End At: %s", asctime(tb) );
    LogBuf[strlen( LogBuf )-1] = '\0';
    print_log( LogBuf );

    sprintf( LogBuf, "Total Time %i", time2-time1 );
    print_log( LogBuf );
    print_log( sep_str );

    fclose( LogFilefd );
    MPI_Finalize();
}

