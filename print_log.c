#include "allvars.h"

void print_log( char *log ) {
    fprintf( LogFilefd, "%s\n", log );
    if ( ThisTask == 0 )
        printf( "%s\n", log );
}
