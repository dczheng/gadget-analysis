#include "stdio.h"
#include "stdlib.h"
#include "string.h"

void main( int argc, char **argv ) {

    if ( argc < 5 )
        printf( "ERROR: invalid parameters!\n" );

    char *fn, *buf, *delim=",", *token;
    int nl_f, nl_o, n_fmt, flag, n_col; 
    size_t n_buf;
    FILE *fd;

    fn =  argv[1];
    nl_o = atoi( argv[argc-2] );
    n_fmt = atoi( argv[argc-1] );

    printf( "File: %s\n"
            "Output lines: %i\n"
            "Output width: %i\n",
            fn, nl_o, n_fmt );

    if ( NULL == ( fd = fopen( fn, "r" )) ) {
        printf( "ERROR: can't open file `%s`\n", fn );
    }

    n_buf = 100;
    buf = malloc( n_buf );

    while( !feof(fd) ){
        getline( &buf, &n_buf, fd );
        token = strtok( buf, delim );
        printf( "%i %s", n_buf, buf );
        token = strtok( buf, delim );
        n_col;
        while( token ) {
            token = strtok( NULL, delim );
        }
        break;
    }
    free( buf );
    fclose( fd );

}
