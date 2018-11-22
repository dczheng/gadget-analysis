#include "allvars.h"

#define MAXTAGS 300
#define REAL 1
#define STRING 2
#define INT 3

#define ADD_PARAR( a ) {\
    strcpy( tag[nt], &((#a)[4]) );\
    addr[nt] = &a;\
    id[nt++] = REAL;\
}
#define ADD_PARAI( a ) {\
    strcpy( tag[nt], &((#a)[4]) );\
    addr[nt] = &a;\
    id[nt++] = INT;\
}

#define ADD_PARAS( a ) {\
    strcpy( tag[nt], &((#a)[4]) );\
    addr[nt] = a;\
    id[nt++] = STRING;\
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
            endrun0( "Faile to Open Parameter file %s\n", fn );
        }

        if ( sizeof( long ) != 8 )
            endruns( "Type `long` is no 64 bit on this platform. Stopping." );

        nt = 0;

        ADD_PARAS( All.FilePrefix );
        ADD_PARAS( All.FoFDir     );
        ADD_PARAS( All.RadDir     );
        ADD_PARAS( All.GroupDir   );

        ADD_PARAR( All.UnitMass_in_g            );
        ADD_PARAR( All.UnitLength_in_cm         );
        ADD_PARAR( All.UnitVelocity_in_cm_per_s );
        ADD_PARAR( All.Freq                     );
        ADD_PARAR( All.SofteningGas             );
        ADD_PARAR( All.SofteningHalo            );
        ADD_PARAR( All.SofteningDisk            );
        ADD_PARAR( All.SofteningBulge           );
        ADD_PARAR( All.SofteningStar            );
        ADD_PARAR( All.SofteningBndry           );
        ADD_PARAR( All.StartX                   );
        ADD_PARAR( All.StartY                   );
        ADD_PARAR( All.StartZ                   );
        ADD_PARAR( All.EndX                     );
        ADD_PARAR( All.EndY                     );
        ADD_PARAR( All.EndZ                     );
        ADD_PARAR( All.TreeAllocFactor          );
        ADD_PARAR( All.LinkLength               );
        ADD_PARAR( All.OmegaBaryon              );
        ADD_PARAR( All.ConvSigma                );
        ADD_PARAR( All.Sigma8                   );
        ADD_PARAR( All.GroupMassMin             );
        ADD_PARAR( All.MFMmin                   );
        ADD_PARAR( All.MFMmax                   );
        ADD_PARAR( All.MFMSplit                 );
        ADD_PARAR( All.NuMin                    );
        ADD_PARAR( All.NuMax                    );
        ADD_PARAR( All.QMin                     );
        ADD_PARAR( All.QMax                     );

        ADD_PARAI( All.NumFiles                 );
        ADD_PARAI( All.PicSize                  );
        ADD_PARAI( All.StartSnapIndex           );
        ADD_PARAI( All.ProjectDirection         );
        ADD_PARAI( All.KernelN                  );
        ADD_PARAI( All.FoFMinLen                );
        ADD_PARAI( All.TreePartType             );
        ADD_PARAI( All.ConvN                    );
        ADD_PARAI( All.FoF                      );
        ADD_PARAI( All.GroupIndexMin            );
        ADD_PARAI( All.GroupIndexMax            );
        ADD_PARAI( All.GasState                 );
        ADD_PARAI( All.GasDensity               );
        ADD_PARAI( All.GasTemperature           );
        ADD_PARAI( All.KernelInterpolation      );
        ADD_PARAI( All.MpcFlag                  );
        ADD_PARAI( All.Group                    );
        ADD_PARAI( All.GroupSfr                 );
        ADD_PARAI( All.GroupTemp                );
        ADD_PARAI( All.GroupB                   );
        ADD_PARAI( All.GroupMach                );
        ADD_PARAI( All.GroupHge                 );
        ADD_PARAI( All.GroupRad                 );
        ADD_PARAI( All.ReadMach                 );
        ADD_PARAI( All.ReadTemp                 );
        ADD_PARAI( All.ReadHge                  );
        ADD_PARAI( All.ReadCr                   );
        ADD_PARAI( All.ReadB                    );
        ADD_PARAI( All.TotSpec                  );
        ADD_PARAI( All.ReadSfr                  );
        ADD_PARAI( All.GroupSpec                );
        ADD_PARAI( All.BPdf                     );
        ADD_PARAI( All.MF                       );
        ADD_PARAI( All.MFBins                   );
        ADD_PARAI( All.RadSpec                  );
        ADD_PARAI( All.NuNum                    );
        ADD_PARAI( All.GroupEleSpec             );
        ADD_PARAI( All.QNum                     );
        ADD_PARAI( All.PowSpec                  );
        ADD_PARAI( All.PowSpecNGrid             );
        ADD_PARAI( All.PowSpecPartType          );
        ADD_PARAI( All.PowSpecBins              );
        ADD_PARAI( All.CrePressurePdf           );


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
            endrun(20181107);
        fclose( fd );
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Bcast( &All, sizeof( struct global_parameters_struct ), MPI_BYTE, 0, MPI_COMM_WORLD );
    put_block_line;
}

