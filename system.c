#include "allvars.h"

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
