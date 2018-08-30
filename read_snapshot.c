#include "allvars.h"
//#define GROUP_ZDEBUG

hid_t hdf5_file, hdf5_group, hdf5_dataset, hdf5_dataspace, hdf5_dataspace_in_file,
       hdf5_dataspace_in_memory, hdf5_type, hdf5_hdf5_type_mem, hdf5_attribute, hdf5_type;
herr_t herr;
hsize_t dims[2], maxdims[2];
int ndims;
void *CommBuffer;

int blockpresent( enum iofields blk, int pt ) {
    switch( blk ) {
        case IO_POS:
        case IO_VEL:
        case IO_ID:
            if ( header.npart[pt] > 0 )
                return 1;
            else
                return 0;
        case IO_MASS:
            if ( (header.npart[pt] > 0 ) && header.mass[pt] == 0 )
                return 1;
            else
                return 0;
        case IO_TEMP:
            if ( All.ReadTemperature == 0 )
                return 0;

        case IO_NE:
        case IO_U:
        case IO_RHO:
            if (( pt == 0 ) && ( header.npart[0] != 0 ))
                return 1;
            else
                return 0;

        case IO_MN:
            if (( pt == 0 ) && ( header.npart[0] != 0 ) && All.MachFlag == 1 )
                return 1;
            else
                return 0;

        case IO_MAG:
        case IO_DIVB:
        case IO_DBDT:
            if (( pt == 0 ) && ( header.npart[0] != 0 ) && All.BFlag == 1)
                return 1;
            else
                return 0;

        case IO_SFR:
            if (( pt == 0 ) && ( header.npart[0] != 0 ) && All.SfrFlag == 1)
                return 1;
            else
                return 0;

        case IO_CRE_C0:
        case IO_CRE_Q0:
        case IO_CRE_E0:
        case IO_CRE_n0:
            if (( pt == 0 ) && ( header.npart[0] != 0 ) && All.HgeFlag == 1)
                return 1;
            else
                return 0;

        case IO_CR_Q0:
        case IO_CR_C0:
        case IO_CR_E0:
        case IO_CR_n0:
            if (( pt == 0 ) && ( header.npart[0] != 0 ) && All.CrFlag == 1)
                return 1;
            else
                return 0;

        case IO_POT:
        case IO_ACCEL:
            return 0;
        default:
            return 0;
    }
}

int get_block_nbytes( enum iofields blk ) {
    int block_nbytes = 0;
    switch ( blk ) {
        case IO_POS:
        case IO_VEL:
        case IO_ACCEL:
        case IO_MAG:
            block_nbytes = 3 * sizeof( OutputFloat );
            break;
        case IO_SFR:
        case IO_DIVB:
        case IO_DBDT:
        case IO_MASS:
        case IO_U:
        case IO_TEMP:
        case IO_RHO:
        case IO_POT:
        case IO_NE:
        case IO_MN:
        case IO_CR_C0:
        case IO_CR_Q0:
        case IO_CR_E0:
        case IO_CR_n0:
        case IO_CRE_C0:
        case IO_CRE_Q0:
        case IO_CRE_E0:
        case IO_CRE_n0:
            block_nbytes = sizeof( OutputFloat );
            break;
        case IO_ID:
            block_nbytes = sizeof( MyIDType  );
            break;
    }
    return block_nbytes;
}

void get_block_dims( int pt, enum iofields blk, hsize_t (*dims)[2] ) {
    switch ( blk ) {
        case IO_POS:
        case IO_VEL:
        case IO_ACCEL:
        case IO_MAG:
            (*dims)[0] = header.npart[pt];
            (*dims)[1] = 3;
            break;
        case IO_SFR:
        case IO_DIVB:
        case IO_DBDT:
        case IO_MASS:
        case IO_U:
        case IO_TEMP:
        case IO_RHO:
        case IO_POT:
        case IO_NE:
        case IO_ID:
        case IO_MN:
        case IO_CR_C0:
        case IO_CR_Q0:
        case IO_CR_E0:
        case IO_CR_n0:
        case IO_CRE_C0:
        case IO_CRE_Q0:
        case IO_CRE_E0:
        case IO_CRE_n0:
            (*dims)[0] = header.npart[pt];
            (*dims)[1] = 1;
            break;
    }
}

void get_dataset_name( enum iofields blk, char *buf ) {
    switch ( blk ) {
        case IO_POS:
            strcpy( buf, "Coordinates" );
            break;
        case IO_VEL:
            strcpy( buf, "Velocities" );
            break;
        case IO_ACCEL:
            strcpy( buf, "Acceleration" );
            break;
        case IO_SFR:
            strcpy( buf, "StarFormationRate" );
            break;
        case IO_MAG:
            strcpy( buf, "MagneticField" );
            break;
        case IO_DIVB:
            strcpy( buf, "DivergenceOfMagneticField" );
            break;
        case IO_DBDT:
            strcpy( buf, "RateOfChangeOfMagneticField" );
            break;
        case IO_MASS:
            strcpy( buf, "Masses" );
            break;
        case IO_U:
            strcpy( buf, "InternalEnergy" );
            break;
        case IO_TEMP:
            strcpy( buf, "Temperature" );
            break;
        case IO_RHO:
            strcpy( buf, "Density" );
            break;
        case IO_POT:
            strcpy( buf, "Potential" );
            break;
        case IO_NE:
            strcpy( buf, "ElectronAbundance" );
            break;
        case IO_ID:
            strcpy( buf, "ParticleIDs" );
            break;
        case IO_MN:
            strcpy( buf, "MachNumber" );
            break;
        case IO_CR_C0:
            strcpy( buf, "CR_C0" );
            break;
        case IO_CR_Q0:
            strcpy( buf, "CR_q0" );
            break;
        case IO_CR_E0:
            strcpy( buf, "CR_E0" );
            break;
        case IO_CR_n0:
            strcpy( buf, "CR_n0" );
            break;
        case IO_CRE_C0:
            strcpy( buf, "CRE_C0" );
            break;
        case IO_CRE_Q0:
            strcpy( buf, "CRE_q0" );
            break;
        case IO_CRE_E0:
            strcpy( buf, "CRE_E0" );
            break;
        case IO_CRE_n0:
            strcpy( buf, "CRE_n0" );
            break;
    }
}

void get_hdf5_native_type( enum iofields blk, hid_t *hdf5_type ) {
    switch ( blk ) {
        case IO_ID:
#ifdef LONGIDS
            *hdf5_type = H5Tcopy(H5T_NATIVE_UINT64);
#else
            *hdf5_type = H5Tcopy(H5T_NATIVE_UINT);
#endif
            break;
        default:
#ifdef OUTPUT_IN_DOUBLEPRECISION
            *hdf5_type = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
            *hdf5_type = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
            break;
    }
}

void empty_buffer( enum iofields blk, int offset, int pt ) {
    int i, j;
    long n;
    OutputFloat *fp;
    MyIDType *ip;
    fp = CommBuffer;
    ip = CommBuffer;
    n = header.npart[pt];
    switch( blk ) {
        case IO_POS:
            for ( i=0; i<n; i++ ){
                for ( j=0; j<3; j++ )
                    P[offset+i].Pos[j] = *fp++;
                P[offset+i].Type = pt;
            }
            break;
        case IO_VEL:
            for ( i=0; i<n; i++ )
                for ( j=0; j<3; j++ )
                    P[offset+i].Vel[j] = *fp++;
            break;
        case IO_ID:
            for ( i=0; i<n; i++ )
                P[offset+i].ID = *ip++;
            break;
        case IO_MASS:
            for ( i=0; i<n; i++ )
                P[offset+i].Mass = *fp++;

            /*
            for ( i=0; i<n; i++ ) {
                if ( pt == 5 )
                    printf( "%g\n", P[offset+i].Mass );
            }
            */

            break;
        case IO_POT:
            for ( i=0; i<n; i++ )
                P[offset+i].Pot = *fp++;
            break;
        case IO_ACCEL:
            for ( i=0; i<n; i++ )
                for ( j=0; j<3; j++ )
                    P[offset+i].Acc[j] = *fp++;
            break;
        case IO_RHO:
            for ( i=0; i<n; i++ )
                SphP[offset+i].Density = *fp++;
            break;
        case IO_SFR:
            for ( i=0; i<n; i++ )
                SphP[offset+i].sfr = *fp++;
            break;
        case IO_MAG:
            for ( i=0; i<n; i++ )
                for ( j=0; j<3; j++ )
                    SphP[offset+i].B[j] = *fp++;
            break;
        case IO_DIVB:
            for ( i=0; i<n; i++ )
                SphP[offset+i].divB = *fp++;
            break;
        case IO_DBDT:
            for ( i=0; i<n; i++ )
                SphP[offset+i].dBdt = *fp++;
            break;
        case IO_U:
            for ( i=0; i<n; i++ )
                SphP[offset+i].u = *fp++;
            break;
        case IO_TEMP:
            for ( i=0; i<n; i++ )
                SphP[offset+i].Temp = *fp++;
            break;
        case IO_NE:
            for ( i=0; i<n; i++ )
                SphP[offset+i].elec = *fp++;
            break;
        case IO_MN:
            for ( i=0; i<n; i++ )
                SphP[offset+i].MachNumber = *fp++;
            break;
        case IO_CR_C0:
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_C0 = *fp++;
            break;
        case IO_CR_Q0:
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_Q0 = *fp++;
            break;
        case IO_CR_E0:
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_E0 = *fp++;
            break;
        case IO_CR_n0:
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_n0 = *fp++;
            break;
        case IO_CRE_C0:
            for ( i=0; i<n; i++ )
                SphP[offset+i].CRE_C0 = *fp++;
            break;
        case IO_CRE_Q0:
            for ( i=0; i<n; i++ )
                SphP[offset+i].CRE_Q0 = *fp++;
            break;
        case IO_CRE_E0:
            for ( i=0; i<n; i++ )
                SphP[offset+i].CRE_E0 = *fp++;
            break;
        case IO_CRE_n0:
            for ( i=0; i<n; i++ )
                SphP[offset+i].CRE_n0 = *fp++;
            break;
        default:
            break;

    }
}

void read_header( char *fn ) {
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

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_IC_Info");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_ic_info);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Sfr");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Cooling");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_StellarAge");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Metals");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Feedback");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_DoublePrecision");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision);
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

    hdf5_attribute = H5Aopen_name(hdf5_group, "Redshift");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
    H5Aclose(hdf5_attribute);

    H5Gclose(hdf5_group);
    H5Fclose(hdf5_file);
}

void show_header( struct io_header header ) {
    int i;

    put_block_line;
    writelog( "header Info:\n" );

    writelog( "%-25s: ", "npart" );
    for ( i=0; i<6; i++ )
        writelog(  "%i ", header.npart[i] );
    writelog( "\n" );

    writelog(  "%-25s: ", "mass" );
    for ( i=0; i<6; i++ )
        writelog(  "%lf ", header.mass[i] );
    writelog( "\n" );

    writelog( "%-25s: ", "npartTotal" );
    for ( i=0; i<6; i++ )
        writelog( "%i ", header.npartTotal[i] );
    writelog( "\n" );

    writelog( "%-25s: ", "npartTotalHighWord" );
    for ( i=0; i<6; i++ )
        writelog( "%i ", header.npartTotalHighWord[i] );
    writelog( "\n" );

    writelog( "%-25s: %lf\n", "readshift",             header.redshift );
    writelog( "%-25s: %lf\n", "HubbleParam",           header.HubbleParam );
    writelog( "%-25s: %i\n", "flag_sfr",               header.flag_sfr );
    writelog( "%-25s: %i\n", "flag_feedback",          header.flag_feedback );
    writelog( "%-25s: %i\n", "num_files",              header.num_files );
    writelog( "%-25s: %lf\n", "BoxSize",               header.BoxSize );
    writelog( "%-25s: %lf\n", "Omega0",                header.Omega0 );
    writelog( "%-25s: %lf\n", "OmegaLambda",           header.OmegaLambda );
    writelog( "%-25s: %i\n", "flag_stellarge",         header.flag_stellarage );
    writelog( "%-25s: %i\n", "flag_metals",            header.flag_metals );
    writelog( "%-25s: %i\n", "flag_entropy_instead_u", header.flag_entropy_instead_u );
    writelog( "%-25s: %i\n", "flag_doubleprecision",   header.flag_doubleprecision );
    writelog( "%-25s: %i\n", "flag_ic_info",           header.flag_ic_info );
    writelog( "%-25s: %f\n", "lpt_scalingfactor",      header.lpt_scalingfactor );
    put_block_line;
}

void write_header( char *fn, struct io_header header ) {

    writelog( "write header To %s\n", fn  );
    hdf5_file = H5Fcreate( fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
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
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_DoublePrecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_IC_Info", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_ic_info );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );
    H5Gclose( hdf5_group );
    H5Fclose( hdf5_file );
}

void allocate_memory() {
    ms.max_mem = ms.mem = ms.nn = 0;
    mymalloc1( P, NumPart * sizeof( struct particle_data ) );
    mymalloc1( SphP, N_Gas * sizeof( struct sph_particle_data ) );
}

void free_particle_memory() {
    writelog( "free memory ...\n" );
    myfree( P );
    myfree( SphP );
    id_to_index++;
    myfree( id_to_index );
    writelog( "free memory ... done.\n" );
    put_block_line;
}

void find_id() {
    int bits;
    long i, num, offset;
    timer_start();
    writelog( "find id ...\n" );

    for ( bits=0; GENERATIONS > (1<<bits); bits++ );
    writelog( "bits = %i\n", bits );
    writelog( "find gas id ...\n" );
    for ( i=0; i<N_Gas; i++ ){
        P[i].ID <<= bits;
        P[i].ID >>= bits;
    }

    num = get_particle_num( 4 );
    offset = get_particle_offset( 4 );
    if ( num != 0 ) {
        writelog( "Star Particle Offset: %li\n", offset );
        writelog( "find star id ...\n" );
        for ( i=0; i<num; i++ ) {
            P[offset+i].ID <<= bits;
            P[offset+i].ID >>= bits;
        }
    }

    num = get_particle_num( 5 );
    offset = get_particle_offset( 5 );
    if ( num != 0 ) {
        writelog( "Black Hole Particle Offset: %li\n", offset );
        writelog( "find Black Hole id ...\n" );
        for ( i=0; i<num; i++ ) {
            P[offset+i].ID <<= bits;
            P[offset+i].ID >>= bits;
        }
    }

    writelog( "find id ... done.\n" );
    timer_end();
    put_block_line;
}

void construct_id_to_index() {
    long idmax, idmin, i, idn, N;
    timer_start();
    writelog( "construct id to index ...\n" );
    idmax = -1;
    idmin = LONG_MAX;

    for ( i=0; i<NumPart; i++ ) {
        idmax = vmax( idmax, P[i].ID );
        idmin = vmin( idmin, P[i].ID, 0 );
    }

    idn = idmax - idmin + 1;
    writelog( "Min ID: %li, Max ID: %li, ID Num: %li\n", idmin, idmax, idn );

    mymalloc3( id_to_index, idn * sizeof( long ), -1 );

    id_to_index--;

    N = N_Gas + get_particle_num( 1 );

    for ( i=0; i<N; i++ ) {
        id_to_index[ P[i].ID ] = i;
    }

    /*
    long offset, num, id;
    offset = get_particle_offset( 4 );
    num = get_particle_num( 4 );
    for ( i=offset; i<offset+num; i++ ) {
        id = P[i].ID;
        if ( P[id_to_index[id]].ID != id )
            printf( "%li %li %li\n", id, id_to_index[id], P[id_to_index[id]].ID );
    }
    */

    writelog( "construct id to index ... done.\n" );
    timer_end();
    put_block_line;
}

void attach_some_Star_BH_to_Gas() {
    long i, num1, offset1, id, index, n, flag,
         offset2, num2;

    writelog( "attach some star and bh to gas...\n" );

    writelog( "npart: " );
    for( i=0; i<6; i++ )
        writelog( "%li ", get_particle_num( i ) );
    writelog( "\n" );
    writelog( "Total: %li\n", NumPart );

    flag = 0;

    offset1 = get_particle_offset( 4 );
    num1 = get_particle_num( 4 );
    if ( num1 != 0 ) {

        flag = 1;
        for ( i=offset1, n=0; i<offset1+num1; i++ ) {
            id = P[i].ID;
            index = id_to_index[ id ];

            if ( index < 0 )
                continue;

            if ( P[index].Type != 0 )
                endrun( "some error appear." );

            SphP[index].Star_Mass += P[i].Mass;
            n++;
            P[i] = P[ offset1+num1-1 ];
            num1--;
            i--;
        }

        writelog( "attach %li star particle to gas...\n", n );

        NumPart -= n;
        header.npartTotal[4] = ( (num1 << 32) >> 32);
        header.npartTotalHighWord[4] = ( num1 >> 32);

    }


    offset2 = get_particle_offset( 5 );
    num2 = get_particle_num( 5 );
    if ( num2 != 0 ) {

        flag = 1;
        for ( i=offset2, n=0; i<offset2+num2; i++ ) {
            id = P[i].ID;
            index = id_to_index[ id ];

            if ( index < 0 )
                continue;

            if ( P[index].Type != 0 )
                endrun( "some error appear." );

            SphP[index].BH_Mass += P[i].Mass;
            n++;
            P[i] = P[ offset2+num2-1 ];
            num2--;
            i--;
        }

        writelog( "attach %li black hole particle to gas...\n", n );

        NumPart -= n;
        header.npartTotal[5] = ( (num2 << 32) >> 32);
        header.npartTotalHighWord[5] = ( num2 >> 32);

        for ( i=offset2, index=offset1+num1; i<offset2+num2; i++, index++ )
            P[index] = P[i];

    }



    if ( flag == 1 ) {
        writelog( "npart: " );
        for( i=0; i<6; i++ )
            writelog( "%li ", get_particle_num( i ) );
        writelog( "\n" );
        writelog( "Total: %li\n", NumPart );
    }
    else {
        writelog( "there is no star or bh particle to attach gas\n" );
    }

    writelog( "attach some star and bh to gas... done.\n" );
    put_block_line;

}

void check_data( int err ) {

    long i, offset, num;

    writelog( "Check data ...\n" );

    offset = get_particle_offset( 4 );
    num = get_particle_num( 4 );
    printf( "type: 4, offset: %li, num: %li\n", offset, num );

    for ( i=offset; i<offset+num; i++ )
        printf( "%g\n", P[i].Mass );

    offset = get_particle_offset( 5 );
    num = get_particle_num( 5 );
    printf( "type: 5, offset: %li, num: %li\n", offset, num );

    /*
    for ( i=offset; i<offset+num; i++ )
        printf( "%g\n", P[i].Mass );
    */

    for ( i=0; i<NumPart; i++ )
        if ( P[i].Mass == 0 ) {
            //printf( "P[%li].Mass = 0 !!!, Type: %i\n", i, P[i].Type );
            endrun();
        }


    writelog( "Check data ... done.\n" );
    put_block_line;

}

void read_snapshot() {
    int pt, blk, nbytes;
    long i, file, offset, num;
    char file_name[MYFILENAME_MAX], buf[200], buf1[200];
    size_t BufferBytes;
    timer_start();
    writelog( "read_data ...\n" );

#ifdef OUTPUT_IN_DOUBLEPRECISION
    writelog( "******Note: output precision is double ! ******\n" );
#else
    writelog( "******Note: output precision is float ! ******\n" );
#endif

    sprintf( file_name, "%s_%03d.%3i.hdf5", All.FilePrefix, All.StartSnapIndex + ThisTask, 0 );
    if ( All.NumFiles < 2 )
        sprintf( file_name, "%s_%03d.hdf5", All.FilePrefix, All.StartSnapIndex + ThisTask );
    N_Gas = NumPart = 0;
    read_header( file_name );
    /*
    for ( i=0; i<6; i++ ) {
        printf( "%.20f\n", header.mass[i] );
    }
    endrun();
    */
    show_header( header );
    for ( i=0; i<6; i++ ){
        NumPart += get_particle_num( i );
    }
    N_Gas = header.npartTotal[0] + ( ( (long)header.npartTotalHighWord[0] ) << 32 );

    All.BoxSize = header.BoxSize;
    All.HalfBoxSize = All.BoxSize / 2;
    All.RedShift = header.redshift;
    All.Omega0 = header.Omega0;
    All.OmegaLambda = header.OmegaLambda;
    All.HubbleParam = header.HubbleParam;
    writelog( "NumPart = %ld, N_Gas = %ld\n", NumPart, N_Gas );

    allocate_memory();

    BufferBytes = header.npart[1] * sizeof( OutputFloat ) * 3;
    mymalloc1( CommBuffer, BufferBytes );

    for ( blk=0; blk<IO_NBLOCKS; blk++ ) {
        for ( pt=0, offset=0; pt<6; pt++ ) {
            for (file=0; file<All.NumFiles; file++) {
                if ( All.NumFiles < 2 )
                    sprintf( file_name, "%s_%03d.hdf5", All.FilePrefix, All.StartSnapIndex + ThisTask );
                else
                    sprintf( file_name, "%s_%03d.%3li.hdf5",
                            All.FilePrefix, All.StartSnapIndex + ThisTask, file );
                    read_header( file_name );

                    if ( blockpresent( blk, pt ) ) {
                        nbytes = get_block_nbytes( blk );
                        if ( nbytes * header.npart[pt] > BufferBytes ){
                            writelog( "CommBuffer IS TOO SMALL.\n" );
                            myfree( CommBuffer );
                            BufferBytes = nbytes * header.npart[pt];
                            mymalloc1( CommBuffer, BufferBytes );
                        }

                        hdf5_file = H5Fopen( file_name, H5F_ACC_RDWR, H5P_DEFAULT );
                        get_dataset_name( blk, buf );
                        get_hdf5_native_type( blk, &hdf5_type );
                        writelog( "[%i] %8li reading %s ...\n", pt, offset, buf );
                        sprintf( buf1, "PartType%i/%s", pt, buf );
                        hdf5_dataset = H5Dopen( hdf5_file, buf1 );
                        herr = H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, CommBuffer );
                        empty_buffer( blk, offset, pt );
                        H5Dclose( hdf5_dataset );
                        H5Tclose( hdf5_type );
                        H5Fclose( hdf5_file );
                }

                offset += header.npart[pt];
            }
        }
    }
    myfree( CommBuffer );

    //check_data( 212121 );

    for ( pt=0, offset=0; pt<6; pt++ ) {
        num = header.npartTotal[pt] + ( ( ( long )header.npartTotalHighWord[pt] ) << 32 );
        if ( num != 0 && header.mass[pt] != 0 ){
            writelog( "set particle `%d` mass from header\n", pt );
            for ( i=offset; i<offset+num; i++ )
                P[i].Mass = header.mass[pt];
        }
        offset += num;
    }

    for ( i=0; i<N_Gas; i++ )
        SphP[i].BH_Mass = SphP[i].Star_Mass = 0;

    writelog( "read data ... done. \n" );
    timer_end();
    put_block_line;
    find_id();
    construct_id_to_index();
    attach_some_Star_BH_to_Gas();

    //endrun();
}

