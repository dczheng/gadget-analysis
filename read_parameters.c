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
        writelog( "read parameter...\n" );
    if ( ThisTask == 0 ) {
        fd = fopen( fn, "r" );
        if ( NULL == fd ) {
            fprintf( stderr, "Faile to Open Parameter file %s\n", fn );
            endrun( 1 );
        }

        if ( sizeof( long ) != 8 ) {
            printf( "Type `long` is no 64 bit on this platform. Stopping. \n" );
            endrun( 20171207 );
        }

        nt = 0;
        strcpy( tag[nt], "FilePrefix" );
        addr[nt] = All.FilePrefix;
        id[nt++] = STRING;

        strcpy( tag[nt], "NumFiles" );
        addr[nt] = &All.NumFiles;
        id[nt++] = INT;

        strcpy( tag[nt], "PicSize" );
        addr[nt] = &All.PicSize;
        id[nt++] = INT;

        strcpy( tag[nt], "UnitMass_in_g" );
        addr[nt] = &All.UnitMass_in_g;
        id[nt++] = REAL;

        strcpy( tag[nt], "UnitLength_in_cm" );
        addr[nt] = &All.UnitLength_in_cm;
        id[nt++] = REAL;

        strcpy( tag[nt], "UnitVelocity_in_cm_per_s" );
        addr[nt] = &All.UnitVelocity_in_cm_per_s;
        id[nt++] = REAL;

        strcpy( tag[nt], "MpcFlag" );
        addr[nt] = &All.MpcFlag;
        id[nt++] = INT;

        strcpy( tag[nt], "GroupFlag" );
        addr[nt] = &All.GroupFlag;
        id[nt++] = INT;

        strcpy( tag[nt], "MachFlag" );
        addr[nt] = &All.MachFlag;
        id[nt++] = INT;

        strcpy( tag[nt], "FoFRead" );
        addr[nt] = &All.FoFRead;
        id[nt++] = INT;

        strcpy( tag[nt], "SofteningGas" );
        addr[nt] = &All.SofteningTable[0];
        id[nt++] = REAL;

        strcpy( tag[nt], "SofteningHalo" );
        addr[nt] = &All.SofteningTable[1];
        id[nt++] = REAL;

        strcpy( tag[nt], "SofteningDisk" );
        addr[nt] = &All.SofteningTable[2];
        id[nt++] = REAL;

        strcpy( tag[nt], "SofteningBulge" );
        addr[nt] = &All.SofteningTable[3];
        id[nt++] = REAL;

        strcpy( tag[nt], "SofteningStar" );
        addr[nt] = &All.SofteningTable[4];
        id[nt++] = REAL;

        strcpy( tag[nt], "SofteningBndry" );
        addr[nt] = &All.SofteningTable[5];
        id[nt++] = REAL;

        strcpy( tag[nt], "Alpha" );
        addr[nt] = &All.Alpha;
        id[nt++] = REAL;

        strcpy( tag[nt], "StartSnapIndex" );
        addr[nt] = &All.StartSnapIndex;
        id[nt++] = INT;

        strcpy( tag[nt], "ProjectDirection" );
        addr[nt] = &All.ProjectDirection;
        id[nt++] = INT;

        strcpy( tag[nt], "KernelN" );
        addr[nt] = &All.KernelN;
        id[nt++] = INT;

        strcpy( tag[nt], "StartX" );
        addr[nt] = &All.Start[0];
        id[nt++] = REAL;

        strcpy( tag[nt], "EndX" );
        addr[nt] = &All.End[0];
        id[nt++] = REAL;

        strcpy( tag[nt], "StartY" );
        addr[nt] = &All.Start[1];
        id[nt++] = REAL;

        strcpy( tag[nt], "EndY" );
        addr[nt] = &All.End[1];
        id[nt++] = REAL;

        strcpy( tag[nt], "StartZ" );
        addr[nt] = &All.Start[2];
        id[nt++] = REAL;

        strcpy( tag[nt], "EndZ" );
        addr[nt] = &All.End[2];
        id[nt++] = REAL;

        strcpy( tag[nt], "HgeFlag" );
        addr[nt] = &All.HgeFlag;
        id[nt++] = INT;

        strcpy( tag[nt], "CrFlag" );
        addr[nt] = &All.CrFlag;
        id[nt++] = INT;

        strcpy( tag[nt], "BFlag" );
        addr[nt] = &All.BFlag;
        id[nt++] = INT;

        strcpy( tag[nt], "GasState" );
        addr[nt] = &All.GasState;
        id[nt++] = INT;

        strcpy( tag[nt], "ReadTemperature" );
        addr[nt] = &All.ReadTemperature;
        id[nt++] = INT;

        strcpy( tag[nt], "GasDensity" );
        addr[nt] = &All.GasDensity;
        id[nt++] = INT;

        strcpy( tag[nt], "GasTemperature" );
        addr[nt] = &All.GasTemperature;
        id[nt++] = INT;

        strcpy( tag[nt], "KernelInterpolation" );
        addr[nt] = &All.KernelInterpolation;
        id[nt++] = INT;

        strcpy( tag[nt], "TreeAllocFactor" );
        addr[nt] = &All.TreeAllocFactor;
        id[nt++] = REAL;

        strcpy( tag[nt], "LinkLength" );
        addr[nt] = &All.LinkLength;
        id[nt++] = REAL;

        strcpy( tag[nt], "FoFMinLen" );
        addr[nt] = &All.FoFMinLen;
        id[nt++] = INT;

        strcpy( tag[nt], "TreePartType" );
        addr[nt] = &All.TreePartType;
        id[nt++] = INT;

        strcpy( tag[nt], "OmegaBaryon" );
        addr[nt] = &All.OmegaBaryon;
        id[nt++] = REAL;

        strcpy( tag[nt], "ConvN" );
        addr[nt] = &All.ConvN;
        id[nt++] = INT;

        strcpy( tag[nt], "ConvFlag" );
        addr[nt] = &All.ConvFlag;
        id[nt++] = INT;

        strcpy( tag[nt], "FoF" );
        addr[nt] = &All.FoF;
        id[nt++] = INT;

        strcpy( tag[nt], "MF" );
        addr[nt] = &All.MF;
        id[nt++] = INT;

        strcpy( tag[nt], "ConvSigma" );
        addr[nt] = &All.ConvSigma;
        id[nt++] = REAL;

        strcpy( tag[nt], "FoFPrefix" );
        addr[nt] = &All.FoFPrefix;
        id[nt++] = STRING;

        strcpy( tag[nt], "GroupIndex" );
        addr[nt] = &All.GroupIndex;
        id[nt++] = INT;

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
                        writelog( "%-35s: %g\n", buf1, *((double*)addr[j]) );
                        break;
                    case INT:
                        *( (int*)addr[j] ) = atoi( buf2 );
                        writelog( "%-35s: %d\n", buf1, *((int*)addr[j]) );
                        break;
                    case STRING:
                        strcpy( (char*)addr[j], buf2 );
                        writelog( "%-35s: %s\n", buf1, buf2 );
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
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Bcast( &All, sizeof( struct global_parameters_struct ), MPI_BYTE, 0, MPI_COMM_WORLD );
    put_block_line;
}

