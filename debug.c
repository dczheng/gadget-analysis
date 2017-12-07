#include "allvars.h"

#ifdef DEBUG
void print_vars() {
    fprintf( stdout, "debug_i = %i\n", debug_i );
    fprintf( stdout, "debug_l = %li\n", debug_l );
    fprintf( stdout, "debug_f = %f\n", debug_f );
    fprintf( stdout, "debug_d = %lf\n", debug_d );
    fprintf( stdout, "debug_s = %s\n", debug_s );
}

void signal_hander( int sig ) {
    print_log( sep_str );
    switch ( sig ) {
        case SIGSEGV:
            fprintf( stderr, "Get a signal: SIGSEGV on process %d.\n", ThisTask);
            print_vars();
            print_log( sep_str );
            endrun( 3 );
            break;
    }
}
#endif
