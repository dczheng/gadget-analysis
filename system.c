#include "allvars.h"

double t0, t1;

void create_dir( char *s ) {
    if ( ThisTask == 0 ){
        if ( access( s, 0 ) == -1 ){
            writelog( "create directory `%s` by task 0\n", s );
            if ( mkdir( s, 0755) == -1 ){
                printf( "failed create directory %s.\n", s );
                endrun( 20180424 );
            }
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
}

double second() {
    return ( (double) clock() / CLOCKS_PER_SEC );
}
void endrun( int ierr ) {
    fprintf( stderr, "EXIT CODE: %i\n", ierr );
    MPI_Abort( MPI_COMM_WORLD, ierr );
    exit( ierr );
}
