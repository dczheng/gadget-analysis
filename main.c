#include "allvars.h"

void endrun( int ierr ) {
    fprintf( stderr, "EXIT CODE: %i\n", ierr );
    exit( ierr );
}

void read_para( char *fn ) {
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
                    printf( "%-35s%g\n", buf1, *((double*)addr[j]) );
                    break;
                case INT:
                    *( (int*)addr[j] ) = atoi( buf2 );
                    printf( "%-35s%d\n", buf1, *((int*)addr[j]) );
                    break;
                case STRING:
                    strcpy( (char*)addr[j], buf2 );
                    printf( "%-35s%s\n", buf1, buf2 );
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

void set_units() {
    fputs( sep_str, stdout );
    fprintf( stdout, "Set Units... \n" );
    UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    UnitDensity_in_cgs = UnitMass_in_g / pow( UnitLength_in_cm, 3 );
    UnitEnergy_in_cgs = UnitMass_in_g * pow( UnitLength_in_cm,2 ) / pow( UnitTime_in_s, 2 );
    UnitTime_in_Megayears = UnitTime_in_s / SEC_PER_MEGAYEAR;
    printf( "UnitMass_in_g = %g\n", UnitMass_in_g );
    printf( "UnitTime_in_s = %g\n", UnitTime_in_s );
    printf( "UnitLength_in_cm = %g\n", UnitLength_in_cm );
    printf( "UnitDensity_in_cgs = %g\n", UnitDensity_in_cgs );
    printf( "UnitEnergy_in_cgs = %g\n", UnitEnergy_in_cgs );
    printf( "UnitVelocity_in_cm_per_s = %g\n", UnitVelocity_in_cm_per_s );
    printf( "UnitTime_in_Megayears = %g\n", UnitTime_in_Megayears );
    fputs( sep_str, stdout );
}

void init_sep_str() {
    memset( sep_str, '-', SEP_LEN-1 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';
}

int main( int argc, char *argv[] ){
    int i;
    time_t time1, time2;
    struct tm *tb;
    char tmp[100];
    if ( argc < 2 ) {
        fprintf( stderr, "Parameter file is required on command line!\n " );
        endrun( 1 );
    }
#ifdef DEBUG
    signal( SIGSEGV, signal_hander );
#endif
    time1 = time( NULL );
    tb = localtime( &time1 );
    fprintf( stdout, "Start At: %s", asctime(tb) );
    init_sep_str();
    read_para( argv[1] );
    set_units();
    read_all_data();
    //group_analysis();
    hg_electrons_analysis();
    free_all_memory();
    time2 = time( NULL );
    tb = localtime( &time2 );
    fprintf( stdout, "End At: %s", asctime(tb) );
    fprintf( stdout, "Total Time %i\n", time2-time1 );
}

