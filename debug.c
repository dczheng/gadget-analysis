#include "allvars.h"
#include "protos.h"
#ifdef ZDEBUG
void signal_hander( int s )
{
    int i, ii, k, kk;
    char buf[1000];
    for ( i=0; i<ZDEBUG_NUM; i++ ) {
        kk = 0;
        for ( k=0; k<=strlen( sig.buf[i] ); k++ )
        {
            buf[kk++] = sig.buf[i][k];
            if ( sig.buf[i][k] == '\n' )
                for ( ii=0; ii<5; ii++ )
                    buf[kk++] = ' ';
        }
        printf( "[%d]: %s\n", i, buf );
    }
    printf( "d: " );
    for ( i=0; i<ZDEBUG_NUM; i++ ) {
        printf( "%g ", sig.d[i] );
    }
    printf( "\nSTOP: %s\n", sig.stop );
    fclose( LogFileFd );
    endrun( 20180507 );
}

void empty_sig_buf() {
    init_sig();
}

void init_sig()
{
    int i;
    sprintf( sig.stop, "******" );
    for ( i=0; i<ZDEBUG_NUM; i++ ) {
        sprintf( sig.buf[i], "******" );
        sig.d[i] = 0;
    }
}

#endif
