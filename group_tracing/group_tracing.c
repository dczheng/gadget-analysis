#include "hdf5.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "unistd.h"

#define SEP_LEN 100

char sep_str[ SEP_LEN ];
char fof_file[ 500 ];
char snapshot_prefix[ 500 ];
char out_file[ 500 ];
char group_id_file[ 500 ];
float *pos, *group_pos;
int group_id[ 1000 ], group_num, *group_len_list, files, Ngroups, tot_len;
long long NumPart;
long long Gas_NumPart;
#ifdef LONGIDS
typedef unsigned long long MyIDType;
#else
typedef unsigned int MyIDType;
#endif

MyIDType *Next, *GroupFirst, *GroupID, *GroupLen, *Id;
hid_t hdf5_file, hdf5_group, hdf5_dataset, hdf5_dataspace_in_file, hdf5_dataspace_in_memory;
hid_t hdf5_dataspace, hdf5_type, hdf5_hdf5_type_mem, hdf5_attribute, herr;
hsize_t dims[2];
int ndims;

struct io_header {
  int npart[6];                        /*!< number of particles of each type in this file */
  double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
                                            stored in the mass-block of the snapshot_prefix file, otherwise they are omitted */
  double time;                         /*!< time of snapshot_prefix file */
  double redshift;                     /*!< redshift of snapshot_prefix file */
  int flag_sfr;                        /*!< flags whether the simulation was including star formation */
  int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];          /*!< total number of particles of each type in this snapshot_prefix. This can be
                                            different from npart if one is dealing with a multi-file snapshot_prefix. */
  int flag_cooling;                    /*!< flags whether cooling was included  */
  int num_files;                       /*!< number of files in multi-file snapshot_prefix */
  double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                       /*!< matter density in units of critical density */
  double OmegaLambda;                  /*!< cosmological constant parameter */
  double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
  int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
  unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  char fill[60];	               /*!< fills to 256 Bytes */
} header;

void read_header( char *fn ) {

    //fprintf( stdout, "Reading header From %s\n", File_Name );
    hid_t hdf5_file, hdf5_group, hdf5_attribute;
    hdf5_file = H5Fopen( fn, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_group = H5Gopen(hdf5_file, "/Header");

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_ThisFile");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, header.npart);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_Total");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_Total_HighWord");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "MassTable");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "BoxSize");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumFilesPerSnapshot");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Cooling");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Feedback");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Metals");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Sfr");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_StellarAge");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "HubbleParam");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Omega0");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "OmegaLambda");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
    H5Aclose(hdf5_attribute);

    H5Gclose(hdf5_group);
    H5Fclose(hdf5_file);
}

void show_header( struct io_header header ) {
    int i;
    fputs( sep_str, stdout );
    fputs( "header Info: \n", stdout );
    fprintf( stdout, "%-25s: ", "npart" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npart[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "mass" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%lf ", header.mass[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "npartTotal" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npartTotal[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "npartTotalHighWord" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npartTotalHighWord[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: %lf\n", "readshift",             header.redshift );
    fprintf( stdout, "%-25s: %lf\n", "time",                  header.time );
    fprintf( stdout, "%-25s: %i\n", "flag_sfr",               header.flag_sfr );
    fprintf( stdout, "%-25s: %i\n", "flag_feedback",          header.flag_feedback );
    fprintf( stdout, "%-25s: %i\n", "num_files",              header.num_files );
    fprintf( stdout, "%-25s: %lf\n", "BoxSize",               header.BoxSize );
    fprintf( stdout, "%-25s: %lf\n", "Omega0",                header.Omega0 );
    fprintf( stdout, "%-25s: %lf\n", "OmegaLambda",           header.OmegaLambda );
    fprintf( stdout, "%-25s: %lf\n", "Hubbleparam",           header.HubbleParam );
    fprintf( stdout, "%-25s: %i\n", "flag_stellarge",         header.flag_stellarage );
    fprintf( stdout, "%-25s: %i\n", "flag_metals",            header.flag_metals );
    fprintf( stdout, "%-25s: %i\n", "flag_entropy_instead_u", header.flag_entropy_instead_u );
    fputs( sep_str, stdout );
}

void endrun( int ierr ) {
    fprintf( stderr, "Exit Code: %i\n", ierr );
    exit( ierr );
}

void read_parameter_file(char *fname) {
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300
    FILE *fd, *fdout;
    char buf[200], buf1[200], buf2[200], buf3[400];
    int i, j, nt;
    int id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];
    int  errorFlag = 0;
    nt = 0;
    strcpy(tag[nt], "FofFile");
    addr[nt] = &fof_file;
    id[nt++] = STRING;
    strcpy(tag[nt], "SnapPrefix");
    addr[nt] = &snapshot_prefix;
    id[nt++] = STRING;
    strcpy(tag[nt], "Files");
    addr[nt] = &files;
    id[nt++] = INT;
    strcpy(tag[nt], "OutputFile");
    addr[nt] = &out_file;
    id[nt++] = STRING;
    strcpy(tag[nt], "GroupIdFile");
    addr[nt] = &group_id_file;
    id[nt++] = STRING;

    if(!(fd = fopen(fname, "r")))
    {
        printf("Parameter file %s not found.\n", fname);
        endrun( 1 );
    }
    sprintf(buf, "%s%s", fname, "-usedvalues");
    if(!(fdout = fopen(buf, "w")))
    {
        printf("error opening file '%s' \n", buf);
        endrun( 2 );
    }
    while(!feof(fd))
    {
        *buf = 0;
        fgets(buf, 200, fd);
        if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
            continue;
        if(buf1[0] == '%')
            continue;
        for(i = 0, j = -1; i < nt; i++)
            if(strcmp(buf1, tag[i]) == 0)
            {
                j = i;
                tag[i][0] = 0;
                break;
            }
        if(j >= 0)
        {
            switch (id[j])
            {
                case DOUBLE:
                    *((double *) addr[j]) = atof(buf2);
                    fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
                    break;
                case STRING:
                    strcpy(addr[j], buf2);
                    fprintf(fdout, "%-35s%s\n", buf1, buf2);
                    break;
                case INT:
                    *((int *) addr[j]) = atoi(buf2);
                    fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
                    break;
            }
        }
        else
        {
            fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
                    fname, buf1);
            errorFlag = 1;
        }
    }
    fclose(fd);
    fclose(fdout);
    if(errorFlag)
        for(i = 0; i < nt; i++)
        {
            if(*tag[i])
            {
                printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
                errorFlag = 1;
            }
        }
#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
}

void read_group_id( char *fname ) {
    FILE *fd;
    char buf[500];
    int i;
    fd = fopen( fname, "r" );
    if ( !fd ) {
        fprintf( stderr, "Failed to open file: %s\n", fname );
        endrun( 3 );
    }
    group_num = 0;
    fgets( buf, 500, fd );
    while ( !feof( fd ) ) {
        sscanf( buf, "%d", group_id + group_num );
        group_num ++;
        fgets( buf, 500, fd );
    }
    printf( "group num: %i\n", group_num );
    for ( i=0; i<group_num; i++ ) {
        printf( "%i  ", group_id[i] );
    }
    printf( "\n" );
    fclose( fd );
}

void loadpositions_hdf5( char *fname, int files ) {
    char fn[200];
    long long i, pc;
    for ( i=0, pc=0; i<files; i++ ){
        if(files>1)
            sprintf(fn,"%s.%d.hdf5",fname,i);
        else
            sprintf( fn, "%s.hdf5", fname);
        read_header( fn );
        if ( i==0 ) {
            NumPart=header.npartTotalHighWord[1];
            Gas_NumPart=header.npartTotalHighWord[0];
            NumPart = ( NumPart<<32 ) + header.npartTotal[1];
            Gas_NumPart = ( Gas_NumPart<<32 ) + header.npartTotal[0];
            pos = ( float* )malloc( sizeof( float ) * NumPart * 3 );
            pos -= 3;
            Id = ( MyIDType* )malloc( sizeof( MyIDType ) * NumPart );
            Id--;
        }
        hdf5_file = H5Fopen( fn, H5F_ACC_RDWR, H5P_DEFAULT );
        hdf5_type = H5Tcopy( H5T_NATIVE_FLOAT );
        hdf5_dataset = H5Dopen( hdf5_file, "PartType1/Coordinates" );
        herr = H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos+3+pc*3  );
        H5Dclose( hdf5_dataset );
        H5Tclose( hdf5_type );
#ifdef LONGIDS
        hdf5_type = H5Tcopy( H5T_NATIVE_UINT64 );
#else
        hdf5_type = H5Tcopy( H5T_NATIVE_UINT );
#endif
        hdf5_dataset = H5Dopen( hdf5_file, "PartType1/ParticleIDs" );
        herr = H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, Id+pc+1  );
        H5Dclose( hdf5_dataset );
        H5Tclose( hdf5_type );
        pc += header.npart[1];
        H5Fclose( hdf5_file );
    }
    show_header( header );
}

void loadfof( char *fname ) {
    Next = ( MyIDType* ) malloc( sizeof( MyIDType ) * NumPart );
    Next --;
    ndims = 1;
    dims[0] = NumPart;
    hdf5_file = H5Fopen( fname, H5F_ACC_RDWR, H5P_DEFAULT );
#ifdef LONGIDS
    hdf5_type = H5Tcopy( H5T_NATIVE_UINT64 );
#else
    hdf5_type = H5Tcopy( H5T_NATIVE_UINT );
#endif
    hdf5_dataset = H5Dopen( hdf5_file, "/Next_List" );
    herr = H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, Next+1 );
    H5Dclose( hdf5_dataset );

    hdf5_attribute = H5Aopen_name(hdf5_file, "Num_Groups_Above_Min_Len");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, &Ngroups );
    H5Aclose(hdf5_attribute);
    printf( "Ngroups: %d\n", Ngroups );

    GroupLen = ( MyIDType* ) malloc( sizeof( MyIDType ) * Ngroups );
    GroupFirst = ( MyIDType* ) malloc( sizeof( MyIDType ) * Ngroups );
    GroupID = ( MyIDType* ) malloc( sizeof( MyIDType ) * Ngroups );
    GroupLen--;
    GroupFirst--;
    GroupID--;
    dims[0] = Ngroups;

    hdf5_dataset = H5Dopen( hdf5_file, "/Group_Above_Min_Len_Len" );
    herr = H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, GroupLen+1 );
    H5Dclose( hdf5_dataset );

    hdf5_dataset = H5Dopen( hdf5_file, "/Group_Above_Min_Len_First" );
    herr = H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, GroupFirst+1 );
    H5Dclose( hdf5_dataset );

    hdf5_dataset = H5Dopen( hdf5_file, "/Group_Above_Min_Len_ID" );
    herr = H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, GroupID+1 );
    H5Dclose( hdf5_dataset );
}

void reordering() {
     MyIDType i,j;
    float xyzsave[3], xyzsource[3];
    MyIDType   idsource, idsave, dest;
    MyIDType  signal;
    MyIDType minid,maxid;
    fputs( sep_str, stdout );
    printf("reordering....\n");
    for(i=1,minid=Id[1],maxid=Id[1]; i<=NumPart; i++) {
        if(Id[i] < minid)
            minid= Id[i];
        if(Id[i]>maxid)
            maxid=Id[i];
        }
    printf("minid= %ld  maxid=%ld  maxid-minid=%ld\n", minid,maxid,maxid-minid);
    if ( Gas_NumPart>0 ) {
        printf( "Contain gas, reset id ...\n" );
        for ( i=1; i<=NumPart; i++ ) {
            Id[i] -= Gas_NumPart;
        }
        minid -= Gas_NumPart;
        maxid -= Gas_NumPart;
        printf("minid= %ld  maxid=%ld  maxid-minid=%ld\n", minid,maxid,maxid-minid);
    }
    if ( minid == 0 ) {
        printf( "minid is 0, reset id ...\n" );
        for ( i=1; i<=NumPart; i++ )
            Id[i] ++;
        minid ++;
        maxid ++;
        printf("minid= %ld  maxid=%ld  maxid-minid=%ld\n", minid,maxid,maxid-minid);
    }
    fflush( stdout );
    for(i=1, signal=0; i<=NumPart; i++)
    {
        //printf("%d\n",i);
        //  printf("----%d\n",(signal/100)*NumPart);
        /*
           if(i>(signal/100)*NumPart)
           {
           printf("%d\n",2);
        //	  printf("----%d\n",signal%10);
        if((signal%10)==0)
        {
        printf("\n%d",signal);
        printf("%d",3);
        }
        else
        printf(".");
        fflush(stdout);
        signal++;
        printf("%d",signal);
        }
        */
        if(Id[i] != i)
        {
            for(j=0;j<3;j++)
              xyzsource[j]=pos[i*3+j];
            idsource=Id[i];
            dest=Id[i];
            do
            {
                for(j=0;j<3;j++)
                  xyzsave[j]=pos[dest*3+j];
                idsave=Id[dest];
                for(j=0;j<3;j++)
                  pos[dest*3+j]= xyzsource[j];
                Id[dest]= idsource;
                if(dest == i)
                  break;
                for(j=0;j<3;j++)
                  xyzsource[j]= xyzsave[j];
                idsource=idsave;
                dest=idsource;
            }
            while(1);
        }
    }
    printf("done.\n");
    Id++;
    free(Id);
    printf("space for particle ID freed\n");
    fputs( sep_str, stdout );
    fflush( stdout );
}

void write_hdf5_header() {

    fprintf(  stdout, "write header To %s\n", out_file  );
    hdf5_file = H5Fcreate( out_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    //hdf5_file = H5Fopen( fn, H5F_ACC_RDWR, H5P_DEFAULT );
    hdf5_group = H5Gcreate( hdf5_file, "/Header", 0 );
    hsize_t adim[1] = { 6 };

    hdf5_dataspace = H5Screate( H5S_SIMPLE );
    H5Sset_extent_simple( hdf5_dataspace, 1, adim, NULL );
    hdf5_attribute = H5Acreate( hdf5_group, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_UINT, header.npart );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SIMPLE );
    H5Sset_extent_simple( hdf5_dataspace, 1, adim, NULL );
    hdf5_attribute = H5Acreate( hdf5_group, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SIMPLE );
    H5Sset_extent_simple( hdf5_dataspace, 1, adim, NULL );
    hdf5_attribute = H5Acreate( hdf5_group, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SIMPLE );
    H5Sset_extent_simple( hdf5_dataspace, 1, adim, NULL );
    hdf5_attribute = H5Acreate( hdf5_group, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.num_files );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0 );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_Entropy_ICs", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_entropy_instead_u );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    H5Gclose( hdf5_group );
    H5Fclose( hdf5_file );

}

void write_hdf5_file() {
    long pt,pc;
    char buf[50];
    long ndims, i, j, k;
    float *pos0, *pos1;
    printf( "write to %s... \n", out_file );
    hdf5_file = H5Fopen( out_file, H5F_ACC_RDWR, H5P_DEFAULT );

    pt = 1;
    sprintf( buf, "/PartType%i", pt );
    hdf5_group = H5Gcreate( hdf5_file, buf, 0 );
    ndims = 2;
    dims[0] = header.npart[pt];
    dims[1] = 3;
    hdf5_type = H5Tcopy( H5T_NATIVE_FLOAT );
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    hdf5_dataset = H5Dcreate( hdf5_group, "Coordinates", hdf5_type,
            hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL,
            H5P_DEFAULT, pos+1 );
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );
    H5Dclose( hdf5_dataset );
    H5Gclose( hdf5_group );

    pt = 5;
    sprintf( buf, "/PartType%i", pt );
    hdf5_group = H5Gcreate( hdf5_file, buf, 0 );
    ndims = 2;
    dims[0] = header.npart[pt];
    dims[1] = 3;
    hdf5_type = H5Tcopy( H5T_NATIVE_FLOAT );
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    hdf5_dataset = H5Dcreate( hdf5_group, "Coordinates", hdf5_type,
            hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL,
            H5P_DEFAULT, group_pos );
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );
    H5Dclose( hdf5_dataset );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Group_Num", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &Ngroups );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    ndims = 1;
    dims[0] = group_num;
    hdf5_type = H5Tcopy( H5T_NATIVE_UINT );
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    hdf5_dataset = H5Dcreate( hdf5_group, "Group_Len_List", hdf5_type,
            hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL,
            H5P_DEFAULT, group_len_list );
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );
    H5Dclose( hdf5_dataset );

    ndims = 1;
    dims[0] = group_num;
    hdf5_type = H5Tcopy( H5T_NATIVE_UINT );
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    hdf5_dataset = H5Dcreate( hdf5_group, "Group_Len_Id", hdf5_type,
            hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL,
            H5P_DEFAULT, group_id );
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );
    H5Dclose( hdf5_dataset );

    H5Gclose( hdf5_group  );


    H5Fclose( hdf5_file );
}

void main( int argc, char** argv ) {
    if ( argc < 2 ) {
        fprintf( stderr, "Too few command line parameters! \n" );
        return;
    }
    int i, j, k, g, index;
    memset( sep_str, '-', SEP_LEN-1 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';
    read_parameter_file( argv[1] );
    fputs( sep_str, stdout );
    printf( "fof file          :   %s\n", fof_file );
    printf( "snapshot_prefix   :   %s\n", snapshot_prefix );
    printf( "files             :   %d\n", files );
    printf( "out file          :   %s\n", out_file );
    printf( "group id file     :   %s\n", group_id_file );
    fputs( sep_str, stdout );
    read_group_id( group_id_file );
    fputs( sep_str, stdout );
    loadpositions_hdf5( snapshot_prefix,  files );
    printf( "Numpart:     %lld\n", NumPart );
    printf( "Gas_Numpart: %lld\n", Gas_NumPart );
    fputs( sep_str, stdout );
    reordering();
    loadfof( fof_file );
    fputs( sep_str, stdout );
//----------------------------------------
    tot_len = 0;
    group_len_list = ( int* ) malloc( sizeof( int ) * group_num );
    for ( i=0; i<group_num; i++ ) {
        g = group_id[i];
        group_len_list[i] = GroupLen[g];
        tot_len += GroupLen[g];
        printf( "%i group len: %li\n", g, group_len_list[i] );
    }
    printf( "total len: %li\n", tot_len );
    group_pos = ( float* ) malloc( sizeof( float ) * tot_len * 3 );
    k = 0;
    for ( i=0; i<group_num; i++ ) {
        g = group_id[i];
        index = GroupFirst[g];
        for ( j=0; j<GroupLen[g]; j++ ) {
            group_pos[ k*3+0 ] = pos[ index*3+0 ];
            group_pos[ k*3+1 ] = pos[ index*3+1 ];
            group_pos[ k*3+2 ] = pos[ index*3+2 ];
            k++;
            index = Next[ index ];
        }

    }
    header.npart[1] = header.npartTotal[1];
    header.mass[5]  = header.mass[1];
    header.npart[5] = tot_len;
    header.npartTotal[5] = tot_len;
    write_hdf5_header();
    write_hdf5_file();
//----------------------------------------
    pos += 3;
    free( pos );
    Next++;
    GroupFirst++;
    GroupLen++;
    GroupID++;
    free( Next );
    free( GroupFirst );
    free( GroupLen );
    free( GroupID );
    free( group_pos );
    free( group_len_list );
}


