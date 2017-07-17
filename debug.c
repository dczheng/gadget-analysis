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
    fputs( sep_str, stderr );
    switch ( sig ) {
        case SIGSEGV:
            fprintf( stderr, "Get a signal: SIGSEGV on process %i\n", this_task );
            print_vars();
            fputs( sep_str, stderr );
            end_run( 3 );
            break;
    }
}
#endif
