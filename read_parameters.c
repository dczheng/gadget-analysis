#include "allvars.h"

#define MAXTAGS 300
#define REAL 1
#define STRING 2
#define INT 3

#define ADD_PARA( a, b, c ) {\
    strcpy( tag[nt], a );\
    addr[nt] = b;\
    id[nt++] = c;\
}

void read_parameters( char *fn ) {

    FILE *fd;
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50], buf[200], buf1[200], buf2[200], buf3[200];
        int id[MAXTAGS], nt, i, j, errflag=0;;
        writelog( "read parameter...\n" );
    if ( ThisTask == 0 ) {

        fd = fopen( fn, "r" );

        if ( NULL == fd ){
            endrun1( "Faile to Open Parameter file %s\n", fn );
        }

        if ( sizeof( long ) != 8 )
            endrun( "Type `long` is no 64 bit on this platform. Stopping." );

        nt = 0;

        ADD_PARA( "FilePrefix",                 All.FilePrefix,                  STRING );
        ADD_PARA( "FoFDir",                     All.FoFDir,                      STRING );
        ADD_PARA( "GroupDir",                   All.GroupDir,                    STRING );

        ADD_PARA( "UnitMass_in_g",              &All.UnitMass_in_g,              REAL );
        ADD_PARA( "UnitLength_in_cm",           &All.UnitLength_in_cm,           REAL );
        ADD_PARA( "UnitVelocity_in_cm_per_s",   &All.UnitVelocity_in_cm_per_s,   REAL );
        ADD_PARA( "Freq",                       &All.Freq,                       REAL );
        ADD_PARA( "SofteningGas",               &All.SofteningTable[0],          REAL );
        ADD_PARA( "SofteningHalo",              &All.SofteningTable[1],          REAL );
        ADD_PARA( "SofteningDisk",              &All.SofteningTable[2],          REAL );
        ADD_PARA( "SofteningBulge",             &All.SofteningTable[3],          REAL );
        ADD_PARA( "SofteningStar",              &All.SofteningTable[4],          REAL );
        ADD_PARA( "SofteningBndry",             &All.SofteningTable[5],          REAL );
        ADD_PARA( "Alpha",                      &All.Alpha,                      REAL );
        ADD_PARA( "StartX",                     &All.Start[0],                   REAL );
        ADD_PARA( "StartY",                     &All.Start[1],                   REAL );
        ADD_PARA( "StartZ",                     &All.Start[2],                   REAL );
        ADD_PARA( "EndX",                       &All.End[0],                     REAL );
        ADD_PARA( "EndY",                       &All.End[1],                     REAL );
        ADD_PARA( "EndZ",                       &All.End[2],                     REAL );
        ADD_PARA( "TreeAllocFactor",            &All.TreeAllocFactor,            REAL );
        ADD_PARA( "LinkLength",                 &All.LinkLength,                 REAL );
        ADD_PARA( "OmegaBaryon",                &All.OmegaBaryon,                REAL );
        ADD_PARA( "ConvSigma",                  &All.ConvSigma,                  REAL );
        ADD_PARA( "NuMin",                      &All.NuMin,                      REAL );
        ADD_PARA( "NuMax",                      &All.NuMax,                      REAL );
        ADD_PARA( "GroupMassMin",               &All.GroupMassMin,               REAL );

        ADD_PARA( "NumFiles",                   &All.NumFiles,                   INT );
        ADD_PARA( "PicSize",                    &All.PicSize,                    INT );
        ADD_PARA( "StartSnapIndex",             &All.StartSnapIndex,             INT  );
        ADD_PARA( "ProjectDirection",           &All.ProjectDirection,           INT  );
        ADD_PARA( "KernelN",                    &All.KernelN,                    INT  );
        ADD_PARA( "FoFRead",                    &All.FoFRead,                    INT  );
        ADD_PARA( "FoFMinLen",                  &All.FoFMinLen,                  INT  );
        ADD_PARA( "TreePartType",               &All.TreePartType,               INT  );
        ADD_PARA( "ConvN",                      &All.ConvN,                      INT  );
        ADD_PARA( "FoF",                        &All.FoF,                        INT  );
        ADD_PARA( "NuNum",                      &All.NuNum,                      INT  );
        ADD_PARA( "GroupIndexMin",              &All.GroupIndexMin,              INT  );
        ADD_PARA( "GroupIndexMax",              &All.GroupIndexMax,              INT  );

        ADD_PARA( "GasState",                   &All.GasState,                   INT  );
        ADD_PARA( "GasDensity",                 &All.GasDensity,                 INT  );
        ADD_PARA( "GasTemperature",             &All.GasTemperature,             INT  );
        ADD_PARA( "ReadTemperature",            &All.ReadTemperature,            INT  );
        ADD_PARA( "KernelInterpolation",        &All.KernelInterpolation,        INT  );
        ADD_PARA( "MpcFlag",                    &All.MpcFlag,                    INT  );
        ADD_PARA( "GroupFlag",                  &All.GroupFlag,                  INT  );
        ADD_PARA( "MachFlag",                   &All.MachFlag,                   INT  );
        ADD_PARA( "TempFlag",                   &All.TempFlag,                   INT  );
        ADD_PARA( "MagFlag",                    &All.MagFlag,                    INT  );
        ADD_PARA( "HgeNumDensFlag",             &All.HgeNumDensFlag,             INT  );
        ADD_PARA( "HgeFlag",                    &All.HgeFlag,                    INT  );
        ADD_PARA( "CrFlag",                     &All.CrFlag,                     INT  );
        ADD_PARA( "BFlag",                      &All.BFlag,                      INT  );
        ADD_PARA( "SpecIndexFlag",              &All.SpecIndexFlag,              INT  );
        ADD_PARA( "TotSpecFlag",                &All.TotSpecFlag,                INT  );
        ADD_PARA( "MFFlag",                     &All.MFFlag,                     INT  );
        ADD_PARA( "RadioFlag",                  &All.RadioFlag,                  INT  );
        ADD_PARA( "SfrFlag",                    &All.SfrFlag,                    INT  );
        ADD_PARA( "SpecFlag",                   &All.SpecFlag,                   INT  );

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
            endrun();
        fclose( fd );
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Bcast( &All, sizeof( struct global_parameters_struct ), MPI_BYTE, 0, MPI_COMM_WORLD );
    put_block_line;
}

