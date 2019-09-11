#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hdf5.h"
struct io_header {
    int npart[6];                        /*!< number of particles of each type in this file */
    double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
                                           stored in the mass-block of the snapshot file, otherwise they are omitted */
    double time;                         /*!< time of snapshot file */
    double redshift;                     /*!< redshift of snapshot file */
    int flag_sfr;                        /*!< flags whether the simulation was including star formation */
    int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
    unsigned int npartTotal[6];          /*!< total number of particles of each type in this snapshot. This can be
                                           different from npart if one is dealing with a multi-file snapshot. */
    int flag_cooling;                    /*!< flags whether cooling was included  */
    int num_files;                       /*!< number of files in multi-file snapshot */
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
int NumPart, Ngas;
struct particle_data
{
    float Pos[3];
    float Vel[3];
    float Mass;
    int Type;
    float Rho, U, Temp, Ne;
} *P;
int *Id;
char path[200], filebasename[200],hdf5_file_name[200];
int files;

hid_t hdf5_file, hdf5_group, hdf5_dataset, hdf5_dataspace_in_file, hdf5_dataspace_in_memory, hdf5_dataspace, hdf5_type, hdf5_hdf5_type_mem, hdf5_attribute;
herr_t herr;
hsize_t dims[2], maxdims[2], npoints, precision;

int allocate_memory(void) {
    printf("allocating memory...\n");
    if(!(P = malloc(NumPart * sizeof(struct particle_data))))
    {
        fprintf(stderr, "failed to allocate memory.\n");
        exit(0);
    }
    P--;				/* start with offset 1 */
    if(!(Id = malloc(NumPart * sizeof(int))))
    {
        fprintf(stderr, "failed to allocate memory.\n");
        exit(0);
    }
    Id--;				/* start with offset 1 */
    printf("allocating memory...done\n");
}

int load_snapshot(char *fname, int files) {
    FILE *fd;
    char buf[200];
    int i, j, k, dummy, ntot_withmasses;
    int t, n, off, pc, pc_new, pc_sph;
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
    for(i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
        if(files > 1)
            sprintf(buf, "%s.%d", fname, i);
        else
            sprintf(buf, "%s", fname);
        if(!(fd = fopen(buf, "r")))
        {
            printf("can't open file `%s`\n", buf);
            exit(0);
        }
        printf("reading `%s' ...\n", buf);
        fflush(stdout);
        fread(&dummy, sizeof(dummy), 1, fd);
        fread(&header, sizeof(header), 1, fd);
        fread(&dummy, sizeof(dummy), 1, fd);
        if(files == 1)
        {
            for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
                NumPart += header.npart[k];
            Ngas = header.npart[0];
        }
        else
        {
            for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
                NumPart += header.npartTotal[k];
            Ngas = header.npartTotal[0];
        }
        for(k = 0, ntot_withmasses = 0; k < 6; k++)
        {
            if(header.mass[k] == 0)
                ntot_withmasses += header.npart[k];
        }
        if(i == 0)
            allocate_memory();
        SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header.npart[k]; n++)
            {
                fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;
        SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header.npart[k]; n++)
            {
                fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;
        SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header.npart[k]; n++)
            {
                fread(&Id[pc_new], sizeof(int), 1, fd);
                pc_new++;
            }
        }
        SKIP;
        if(ntot_withmasses > 0)
            SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header.npart[k]; n++)
            {
                P[pc_new].Type = k;
                if(header.mass[k] == 0)
                    fread(&P[pc_new].Mass, sizeof(float), 1, fd);
                else
                    P[pc_new].Mass = header.mass[k];
                pc_new++;
            }
        }
        if(ntot_withmasses > 0)
            SKIP;
        if(header.npart[0] > 0)
        {
            SKIP;
            for(n = 0, pc_sph = pc; n < header.npart[0]; n++)
            {
                fread(&P[pc_sph].U, sizeof(float), 1, fd);
                pc_sph++;
            }
            SKIP;
            SKIP;
            for(n = 0, pc_sph = pc; n < header.npart[0]; n++)
            {
                fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
                pc_sph++;
            }
            SKIP;
            if(header.flag_cooling)
            {
                SKIP;
                for(n = 0, pc_sph = pc; n < header.npart[0]; n++)
                {
                    fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
                    pc_sph++;
                }
                SKIP;
            }
            else
                for(n = 0, pc_sph = pc; n < header.npart[0]; n++)
                {
                    P[pc_sph].Ne = 1.0;
                    pc_sph++;
                }
        }
        fclose(fd);
    }
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
    strcpy(tag[nt], "Path");
    addr[nt] = &path;
    id[nt++] = STRING;
    strcpy(tag[nt], "FileBasename");
    addr[nt] = &filebasename;
    id[nt++] = STRING;
    strcpy(tag[nt], "Files");
    addr[nt] = &files;
    id[nt++] = INT;
    strcpy(tag[nt], "Hdf5_File");
    addr[nt] = &hdf5_file_name;
    id[nt++] = STRING;

    if((fd = fopen(fname, "r")))
    {
        sprintf(buf, "%s%s", fname, "-usedvalues");
        if(!(fdout = fopen(buf, "w")))
        {
            printf("error opening file '%s' \n", buf);
            errorFlag = 1;
        }
        else
        {
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
                            fprintf(stdout, "%-35s%g\n", buf1, *((double *) addr[j]));
                            break;
                        case STRING:
                            strcpy(addr[j], buf2);
                            fprintf(fdout, "%-35s%s\n", buf1, buf2);
                            fprintf(stdout, "%-35s%s\n", buf1, buf2);
                            break;
                        case INT:
                            *((int *) addr[j]) = atoi(buf2);
                            fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
                            fprintf(stdout, "%-35s%d\n", buf1, *((int *) addr[j]));
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
        }
    }
    else
    {
        printf("\nParameter file %s not found.\n\n", fname);
        errorFlag = 2;
    }
    if(errorFlag != 2)
        for(i = 0; i < nt; i++)
        {
            if(*tag[i])
            {
                printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
                errorFlag = 1;
            }
        }
    if(errorFlag)
    {
        exit(0);
    }
#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
}

void show_header( struct io_header header ) {
    int i;
    fputs( "Header2 Info: \n", stdout );
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
    fprintf( stdout, "%-25s: %i\n", "flag_stellarge",         header.flag_stellarage );
    fprintf( stdout, "%-25s: %i\n", "flag_metals",            header.flag_metals );
    fprintf( stdout, "%-25s: %i\n", "flag_entropy_instead_u", header.flag_entropy_instead_u );
}

int reordering(void) {
    int i, j;
    int idsource, idsave, dest, min_id;
    struct particle_data psave, psource;
    printf("reordering....\n");
    min_id = 1000;
    for ( i=1; i<=NumPart; i++ )
        if ( Id[i] < min_id ) {
            min_id = Id[i];
        }
    printf( "minmal id: %d\n", min_id );
    for(i = 1; i <= NumPart; i++)
    {
        if(Id[i] != i)
        {
            psource = P[i];
            idsource = Id[i];
            dest = Id[i];
            do
            {
                psave = P[dest];
                idsave = Id[dest];
                P[dest] = psource;
                Id[dest] = idsource;
                if(dest == i)
                    break;
                psource = psave;
                idsource = idsave;
                dest = idsource;
            }
            while(1);
        }
    }
    printf("done.\n");
    Id++;
    free(Id);
    printf("space for particle ID freed\n");
}

void write_hdf5_header() {

    fprintf(  stdout, "write header To %s\n", hdf5_file_name  );
    hdf5_file = H5Fcreate( hdf5_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
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
    printf( "write to %s... \n", hdf5_file_name );
    hdf5_file = H5Fopen( hdf5_file_name, H5F_ACC_RDWR, H5P_DEFAULT );
    if ( header.npart[0]  != 0 )
        pos0 = ( float* ) malloc( sizeof( float ) * header.npart[0]*3 );
    pos1 = ( float* ) malloc( sizeof( float ) * header.npart[1]*3 );
    /*
       for ( i=0; i<header.npart[1]; i++ ) {
       pos[i*3+0] = P[i+1+header.npart[0]].Pos[0];
       pos[i*3+1] = P[i+1+header.npart[0]].Pos[1];
       pos[i*3+2] = P[i+1+header.npart[0]].Pos[2];
       }
       */
    j = 0;
    k = 0;
    for ( i=1; i<=NumPart; i++ ) {
        if ( Id[i] > header.npart[0] ) {
            pos1[j*3+0] = P[i].Pos[0];
            pos1[j*3+1] = P[i].Pos[1];
            pos1[j*3+2] = P[i].Pos[2];
            j++;
        }
        else {
            if ( header.npart[0]  != 0 )
                pos0[k*3+0] = P[i].Pos[0];
            pos0[k*3+1] = P[i].Pos[1];
            pos0[k*3+2] = P[i].Pos[2];
            k++;
        }
    }
    if ( header.npart[0]  != 0 ) {
        pt = 0;
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
                H5P_DEFAULT, pos0 );
        H5Dclose( hdf5_dataset );
        H5Tclose( hdf5_type );
        H5Sclose( hdf5_dataspace );
        H5Gclose( hdf5_group );
    }

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
            H5P_DEFAULT, pos1 );
    H5Dclose( hdf5_dataset );
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );
    H5Gclose( hdf5_group );

    H5Fclose( hdf5_file );
    if ( header.npart[0]  != 0 )
        free( pos0 );
    free( pos1 );
    /*
       for ( pt=0; pt<2; pt++ ) {
       for ( i=0, pc=0; i<pt; i++ )
       pc += header.npart[i];
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
       H5P_DEFAULT, P[pt].pos );
       H5Tclose( hdf5_type );
       H5Sclose( hdf5_dataspace );

       ndims = 2;
       dims[0] = header.npart[pt];
       dims[1] = 3;
       hdf5_type = H5Tcopy( H5T_NATIVE_FLOAT );
       hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
       hdf5_dataset = H5Dcreate( hdf5_group, "Velocities", hdf5_type,
       hdf5_dataspace, H5P_DEFAULT );
       H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL,
       H5P_DEFAULT, Particle[pt].vel );
       H5Tclose( hdf5_type );
       H5Sclose( hdf5_dataspace );

       ndims = 1;
       dims[0] = header.npart[pt];
       hdf5_type = H5Tcopy( H5T_NATIVE_ULONG );
       hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
       hdf5_dataset = H5Dcreate( hdf5_group, "ParticleIDs", hdf5_type,
       hdf5_dataspace, H5P_DEFAULT );
       H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL,
       H5P_DEFAULT, Particle[pt].id );
       H5Tclose( hdf5_type );
       H5Sclose( hdf5_dataspace );

       H5Gclose( hdf5_group );
       }
       H5Fclose( hdf5_file );
       */
}

int main(int argc, char **argv) {
    if ( argc < 2 ) {
        printf( "Please give a parameter file\n" );
        exit( 0 );
    }
    char input_fname[200];
    int i;
    read_parameter_file( argv[1] );
    sprintf(input_fname, "%s/%s", path, filebasename);
    load_snapshot(input_fname, files);
    printf( "NumPart = %d\n", NumPart );
    //reordering();
    for ( i=0; i<6; i++ ) {
        header.npart[i] = header.npartTotal[i];
    }
    header.num_files = 1;
    show_header( header );
    write_hdf5_header();
    write_hdf5_file();
    P++;
    free( P );
    Id++;
    free( Id );
}
