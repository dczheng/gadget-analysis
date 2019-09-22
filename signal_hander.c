#include "allvars.h"
#include "protos.h"
#ifdef ZDEBUG
#include "string.h"
void signal_hander( int sig )
{
    int i, ii, k, kk, l;
    char buf[ZDEBUG_BUF_LEN];
    FILE *fd;

   sprintf( buf, "signal_hander_info_%i.txt", ThisTask );
   fd = fopen( buf, "w" );

   fprintf( fd, "******************Task %i*******************\n", ThisTask );

   for ( i=0; i<ZDEBUG_NUM; i++ ) {

       kk = 0;
       l = strlen( zsig.buf[i] );

       for ( k=0; k<=l; k++ )
       {

           buf[kk++] = zsig.buf[i][k];

           if ( kk==ZDEBUG_BUF_LEN ) {
               printf( "buffer is too small\n" );
               endrun( 20180916 );
           }

           if ( zsig.buf[i][k] == '\n' )
               for ( ii=0; ii<5; ii++ ) {

                   buf[kk++] = ' ';

                   if ( kk==ZDEBUG_BUF_LEN ) {
                       printf( "buffer is too small\n" );
                       endrun( 20180916 );
                   }

               }

       }

        fprintf( fd, "[%d]: %s\n", i, buf );

   }

   fprintf( fd, "STACK: " );

   for ( i=0; i<zsig.stack_pos; i++ ) {

       if ( (i%5 == 0) && ( i!=0 ) )
           fprintf( fd, "\n       " );

       fprintf( fd, "%g ", zsig.stack[i] );

   }

   fprintf( fd, "\n" );

   fprintf( fd, "STOP: %s\n", zsig.stop );

   fclose( fd );

  //  restart( 0 );
    endrun( 21212121 );

}

void empty_buf_zsig() {

    int i;

    sprintf( zsig.stop, "******" );
    for ( i=0; i<ZDEBUG_NUM; i++ )
        sprintf( zsig.buf[i], "******" );

}

void stack_push_zsig( double a ) {

    zsig.stack[zsig.stack_pos] = a;
    zsig.stack_pos ++;

    if ( zsig.stack_pos == zsig.stack_max ) {

        zsig.stack_max += 1024;
        zsig.stack = realloc( zsig.stack, zsig.stack_max * sizeof( double ) );

    }

}

void stack_pop_zsig() {

    if ( zsig.stack_pos > 0 )
        zsig.stack_pos --;

}

void stack_reset_zsig() {

    zsig.stack_max = 1024;
    zsig.stack_pos = 0;
    free( zsig.stack );
    zsig.stack = malloc( sizeof( double ) * zsig.stack_max );

}

void init_zsig() {

    int i;

    sprintf( zsig.stop, "******" );

    for ( i=0; i<ZDEBUG_NUM; i++ )
        sprintf( zsig.buf[i], "******" );

    zsig.stack_max = 1024;
    zsig.stack_pos = 0;
    zsig.stack = malloc( sizeof( double ) * zsig.stack_max );

}

void free_zsig() {

    free( zsig.stack );

}

MyIDType get_id( int i ) {
    int bits;
    if ( P[i].Type == 0 ) {
        for ( bits=0; GENERATIONS > ( 1<<bits ); bits++ );
        return ( (P[i].ID << bits) >> bits );
    }
    return P[i].ID;
}


#endif
