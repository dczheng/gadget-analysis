#include "allvars.h"

#ifdef DEBUG
void print_vars() {
    int i;
    printf( "debug_l:\n" );
    for ( i=0; i<DEBUG_ARR_LEN; i++ ){
        printf( "%li ", debug_l[i] );
    }
    printf( "\n" );
    printf( "debug_d:\n" );
    for ( i=0; i<DEBUG_ARR_LEN; i++ ){
        printf( "%g ", debug_d[i] );
    }
    printf( "\n" );
    printf( "debug_s\n %s\n", debug_s );
}

void signal_hander( int sig ) {
    writelog( sep_str );
    switch ( sig ) {
        case SIGSEGV:
            fprintf( stderr, "Get a signal: SIGSEGV on process %d.\n", ThisTask);
            print_vars();
            writelog( sep_str );
            endrun( 3 );
            break;
    }
}
#endif
