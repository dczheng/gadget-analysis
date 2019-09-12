#include "allvars.h"

double t0, t1;

void create_dir0( char *s ) {

    if ( ThisTask_Local != 0 )
        return;

    if ( ThisTask_Master == 0 ){
        if ( access( s, 0 ) == -1 ){
            writelog( "`%s` is created by Master 0\n", s );

            if ( mkdir( s, 0755) == -1 )
                endrun0( "failed create directory %s.\n", s );
        }

    }

    MPI_Barrier( MpiComm_Master );

}

void do_sync0( char *s, MPI_Comm comm ) {

    char buf[100];

    if ( NTask == 1 )
        return;

    sprintf( buf, "* synchronization `%s` *", s );
    writelog( "%s\n", buf );
    /*
    int i, n, nn;
    n = strlen( buf );
    nn = strlen( s );
    //printf( "%i\n", nn );

    if ( nn > 0 ) {

        for( i=0; i<n; i++ )
            writelog( "*" );
        writelog( "\n" );
        writelog( "%s\n", buf );
        for( i=0; i<n; i++ )
            writelog( "*" );
        writelog( "\n" );
    }
    */

    MPI_Barrier( comm );
    //endrun( 20181112 );

}

char buf[100];
void do_sync( char *s ) {
    sprintf( buf, "%s [Global]", s );
    do_sync0( buf, MPI_COMM_WORLD );
}

void do_sync_local( char *s ) {
    sprintf( buf, "%s [Local]", s );
    do_sync0( buf, MpiComm_Local );
}

void do_sync_master( char *s ) {
    sprintf( buf, "%s [Master]", s );
    do_sync0( buf, MpiComm_Master );
}


double second() {
    return ( (double) clock() / CLOCKS_PER_SEC );
}

void task_sync_test( char *s ) {

    printf( "[sync test] %s: %i\n", s,ThisTask );
    MPI_Barrier( MPI_COMM_WORLD );

}
