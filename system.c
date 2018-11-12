#include "allvars.h"

double t0, t1;

void create_dir( char *s ) {

    if ( ThisTask == 0 ){
        if ( access( s, 0 ) == -1 ){
            writelog( "create directory `%s` by task 0\n", s );

            if ( mkdir( s, 0755) == -1 )
                endrun0( "failed create directory %s.\n", s );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );

}

void do_sync( char *s ) {

    char buf[100];
    int i, n;

    if ( NTask == 1 )
        return;

    sprintf( buf, "* synchronization after `%s` *", s );

    n = strlen( buf );

    for( i=0; i<n; i++ )
        writelog( "*" );
    writelog( "\n" );
    writelog( "%s\n", buf );
    for( i=0; i<n; i++ )
        writelog( "*" );
    writelog( "\n" );

    MPI_Barrier( MPI_COMM_WORLD );
    //endrun( 20181112 );

}

double second() {
    return ( (double) clock() / CLOCKS_PER_SEC );
}

void task_sync_test() {

    sleep( 2 );
    MPI_Barrier( MPI_COMM_WORLD );
    printf( "%i\n", ThisTask );

}
