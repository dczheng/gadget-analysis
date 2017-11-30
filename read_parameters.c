#include "allvars.h"

void read_parameters( char *fn ) {
#define MAXTAGS 300
#define REAL 1
#define STRING 2
#define INT 3
    FILE *fd;
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50], buf[200], buf1[200], buf2[200], buf3[200];
    int id[MAXTAGS], nt, i, j, errflag=0;;
    fputs( sep_str, stdout );
    fprintf( stdout, "Read Parameter... \n" );
    fd = fopen( fn, "r" );
    if ( NULL == fd ) {
        fprintf( stderr, "Faile to Open Parameter file %s\n", fn );
        endrun( 1 );
    }

    nt = 0;
    strcpy( tag[nt], "FilePrefix" );
    addr[nt] = FilePrefix;
    id[nt++] = STRING;

    strcpy( tag[nt], "GroupDir" );
    addr[nt] = GroupDir;
    id[nt++] = STRING;

    strcpy( tag[nt], "NumFiles" );
    addr[nt] = &NumFiles;
    id[nt++] = INT;

    strcpy( tag[nt], "PicSize" );
    addr[nt] = &PicSize;
    id[nt++] = INT;

    strcpy( tag[nt], "UnitMass_in_g" );
    addr[nt] = &UnitMass_in_g;
    id[nt++] = REAL;

    strcpy( tag[nt], "UnitLength_in_cm" );
    addr[nt] = &UnitLength_in_cm;
    id[nt++] = REAL;

    strcpy( tag[nt], "UnitVelocity_in_cm_per_s" );
    addr[nt] = &UnitVelocity_in_cm_per_s;
    id[nt++] = REAL;

    strcpy( tag[nt], "BufferSize" );
    addr[nt] = &BufferSize;
    id[nt++] = INT;

    strcpy( tag[nt], "MpcFlag" );
    addr[nt] = &MpcFlag;
    id[nt++] = INT;

    strcpy( tag[nt], "SofteningGas" );
    addr[nt] = &SofteningTable[0];
    id[nt++] = REAL;

    strcpy( tag[nt], "SofteningHalo" );
    addr[nt] = &SofteningTable[1];
    id[nt++] = REAL;

    strcpy( tag[nt], "SofteningDisk" );
    addr[nt] = &SofteningTable[2];
    id[nt++] = REAL;

    strcpy( tag[nt], "SofteningBulge" );
    addr[nt] = &SofteningTable[3];
    id[nt++] = REAL;

    strcpy( tag[nt], "SofteningStar" );
    addr[nt] = &SofteningTable[4];
    id[nt++] = REAL;

    strcpy( tag[nt], "SofteningBndry" );
    addr[nt] = &SofteningTable[5];
    id[nt++] = REAL;

    strcpy( tag[nt], "Alpha" );
    addr[nt] = &Alpha;
    id[nt++] = REAL;

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
                    printf( "%-35s: %g\n", buf1, *((double*)addr[j]) );
                    break;
                case INT:
                    *( (int*)addr[j] ) = atoi( buf2 );
                    printf( "%-35s: %d\n", buf1, *((int*)addr[j]) );
                    break;
                case STRING:
                    strcpy( (char*)addr[j], buf2 );
                    printf( "%-35s: %s\n", buf1, buf2 );
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
        endrun( 0 );
    fclose( fd );
}

