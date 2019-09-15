#include "allvars.h"

#define STRING 1
#define INT 2
#define REAL 3

#include "add_params.h"

void read_parameters( char *fn ) {

    FILE *fd, *fd2;
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50], *bname, buf[200], buf1[200], buf2[200], buf3[200];
    int id[MAXTAGS], nt, i, j, errflag=0;

    writelog( "read parameter...\n" );

    if ( ThisTask == 0 ) {

        bname = basename( fn );
        fd2 = myfopen( "w", "./gadget-analysis.log/%s-usedvalues", bname );
        fd = myfopen( "r", fn );

        if ( NULL == fd2 ){
            endrun0( "Faile to Open %s\n", buf );
        }

        if ( NULL == fd ){
            endrun0( "Faile to Open Parameter file %s\n", fn );
        }

        if ( sizeof( long ) != 8 )
            endruns( "Type `long` is no 64 bit on this platform. Stopping." );

        nt = 0;

        ADD_PARAMS();

/*
        for( i=0; i<nt; i++ )
            printf( "%s\n", tag[i] );
*/

        while( !feof( fd ) ) {
            *buf = 0;
            fgets( buf, 200, fd );
            if ( sscanf( buf, "%s%s%s", buf1, buf2, buf3 ) < 2 )
                continue;
            if ( buf1[0] == '%' )
                continue;
            for ( i=0, j=-1; i<nt; i++ )
                if ( strcmp( buf1, tag[i] ) == 0 ) {
                    j = i;
                    tag[i][0] = 0;
                    break;
                }
            if ( j>=0 ) {
                switch ( id[j] ) {
                    case REAL:
                        *( (double*)addr[j] ) = atof( buf2 );
                        writelog( "[D] %-35s: %g\n", buf1, *((double*)addr[j]) );
                        fprintf( fd2, "[D] %-35s %g\n", buf1, *((double*)addr[j]) );
                        break;
                    case INT:
                        *( (int*)addr[j] ) = atoi( buf2 );
                        writelog( "[I] %-35s: %d\n", buf1, *((int*)addr[j]) );
                        fprintf( fd2, "[I] %-35s %d\n", buf1, *((int*)addr[j]) );
                        break;
                    case STRING:
                        strcpy( (char*)addr[j], buf2 );
                        writelog( "[S] %-35s: %s\n", buf1, buf2 );
                        fprintf( fd2, "[S] %-35s %s\n", buf1, buf2 );
                        break;
                }
            }
            else {
                printf( "Error in file %s:  Tag: '%s', not allowed or multiple define.\n", fn, buf1 );
                errflag = 1;
            }
        }
        for ( i=0; i<nt; i++ ) {
            if ( *tag[i] ) {
                printf( "Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fn );
                errflag = 1;
            }
        }
        if ( errflag )
            endrun(20181107);
        fclose( fd );
        fclose( fd2 );
    }


    MPI_Bcast( &All, sizeof( GlobalParams ), MPI_BYTE, 0, MPI_COMM_WORLD );
//    check_flag();
}

