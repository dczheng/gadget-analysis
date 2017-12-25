#include "allvars.h"
//#define GROUP_DEBUG

hid_t hdf5_file, hdf5_group, hdf5_dataset, hdf5_dataspace, hdf5_dataspace_in_file,
       hdf5_dataspace_in_memory, hdf5_type, hdf5_hdf5_type_mem, hdf5_attribute, hdf5_type;
herr_t herr;
hsize_t dims[2], maxdims[2];
int ndims;

int blockpresent( enum iofields blk, int pt ) {
    switch( blk ) {
        case IO_POS:
        case IO_VEL:
        case IO_ID:
            if ( header.npart[pt] > 0 )
                return 1;
        case IO_MASS:
            if ( (header.npart[pt] > 0 ) && header.mass[pt] == 0 )
                return 1;
            else
                return 0;
        case IO_MN:
        case IO_RHO:
            if (( pt == 0 ) && ( header.npart[0] != 0 ))
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
        case IO_U:
        case IO_NE:
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
            block_nbytes = 3 * sizeof( float );
            break;
        case IO_DIVB:
        case IO_DBDT:
        case IO_MASS:
        case IO_U:
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
            block_nbytes = sizeof( float );
            break;
        case IO_ID:
            block_nbytes = sizeof( MyIDType );
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
        case IO_DIVB:
        case IO_DBDT:
        case IO_MASS:
        case IO_U:
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
    int i, j, t;
    long n;
    MyFloat *fp;
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
                SphP[offset+i].Entropy = *fp++;
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
    print_log( sep_str );
    print_log( "header Info:" );

    sprintf( LogBuf, "%-25s: ", "npart" );
    for ( i=0; i<6; i++ )
        sprintf( LogBuf, "%s%li ", LogBuf, header.npart[i] );
    print_log( LogBuf );

    sprintf( LogBuf, "%-25s: ", "mass" );
    for ( i=0; i<6; i++ )
        sprintf( LogBuf, "%s%lf ", LogBuf, header.mass[i] );
    print_log( LogBuf );

    sprintf( LogBuf, "%-25s: ", "npartTotal" );
    for ( i=0; i<6; i++ )
        sprintf( LogBuf, "%s%i ", LogBuf, header.npartTotal[i] );
    print_log( LogBuf );

    sprintf( LogBuf, "%-25s: ", "npartTotalHighWord" );
    for ( i=0; i<6; i++ )
        sprintf( LogBuf, "%s%i ", LogBuf, header.npartTotalHighWord[i] );
    print_log( LogBuf );

    sprintf( LogBuf, "%-25s: %lf", "readshift",             header.redshift );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %lf", "HubbleParam",           header.HubbleParam );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %i", "flag_sfr",               header.flag_sfr );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %i", "flag_feedback",          header.flag_feedback );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %i", "num_files",              header.num_files );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %lf", "BoxSize",               header.BoxSize );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %lf", "Omega0",                header.Omega0 );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %lf", "OmegaLambda",           header.OmegaLambda );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %i", "flag_stellarge",         header.flag_stellarage );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %i", "flag_metals",            header.flag_metals );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %i", "flag_entropy_instead_u", header.flag_entropy_instead_u );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %i", "flag_doubleprecision",   header.flag_doubleprecision );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %i", "flag_ic_info",           header.flag_ic_info );
    print_log( LogBuf );
    sprintf( LogBuf, "%-25s: %f", "lpt_scalingfactor",      header.lpt_scalingfactor );
    print_log( LogBuf );
    print_log( sep_str );
}

void write_header( char *fn, struct io_header header ) {

    sprintf(  LogBuf, "write header To %s", fn  );
    print_log( LogBuf );
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
    double bytes_tot = 0;
    size_t bytes;
    if ( !( P = malloc( bytes = NumPart * sizeof( struct particle_data ) ) ) ) {
        printf( "failed to allocate memory for P ( %g Mb ).\n", bytes / 1024.0 / 1024.0 );
        endrun( 0 );
    }
    bytes_tot += bytes;
    sprintf( LogBuf, "Allocated %g MB for P.", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !( SphP = malloc( bytes = N_Gas * sizeof( struct sph_particle_data ) ) ) ) {
        printf( "failed to allocate memory for SphP ( %g Mb ).\n", bytes / 1024.0 / 1024.0 );
        endrun( 0 );
    }
    bytes_tot += bytes;
    sprintf( LogBuf, "Allocated %g MB for SphP.", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !( CommBuffer = malloc( bytes = All.BufferSize * 1024 * 1024  ) ) ) {
        printf( "failed to allocate memory for CommBuffer ( %g Mb ).\n", bytes / 1024.0 / 1024.0 );
        endrun( 0 );
    }
    bytes_tot += bytes;
    sprintf( LogBuf, "Allocated %g MB for CommBuffer.", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    sprintf( LogBuf, "Allocate total %g MB.", bytes_tot / 1024.0 / 1024.0 );
    print_log( LogBuf );
}

void free_memory() {
    print_log( "free memory ..." );
    free( P );
    free( SphP );
    free( CommBuffer );
    print_log( "free memory ... done." );
    print_log( sep_str );
}

void find_id() {
    int bits;
    long i, num, offset;
    print_log( "find id ..." );

    for ( bits=0; GENERATIONS > (1<<bits); bits++ );
    sprintf( LogBuf, "bits = %i", bits );
    print_log( LogBuf );
    print_log( "find gas id ..." );
    for ( i=0; i<N_Gas; i++ ){
        P[i].ID <<= bits;
        P[i].ID >>= bits;
    }
    num = get_particle_num( 4 );
    offset = get_particle_offset( 4 );
    sprintf( LogBuf, "Start Particle offset %li", offset );
    print_log( LogBuf );
    if ( num != 0 ) {
        print_log( "find gas id ..." );
        for ( i=0; i<num; i++ ) {
            P[offset+i].ID <<= bits;
            P[offset+i].ID >>= bits;
        }
    }
    print_log( "find id ... done." );
    print_log( sep_str );
}

void read_snapshot() {
    int pt, blk, nbytes, rank;
    long i, file, pc, offset, num;
    char file_name[FILENAME_MAX], buf[200], buf1[200];
    print_log( "read_data ..." );
    sprintf( file_name, "%s_%03d.%3i.hdf5", All.FilePrefix, All.StartSnapIndex + ThisTask, 0 );
    if ( All.NumFiles < 2 )
        sprintf( file_name, "%s_%03d.hdf5", All.FilePrefix, All.StartSnapIndex + ThisTask );
    N_Gas = NumPart = 0;
    read_header( file_name );
    for ( i=0; i<6; i++ ){
        NumPart += get_particle_num( i );
    }
    N_Gas = header.npartTotal[0] + ( ( (long long)header.npartTotalHighWord[0] ) << 32 );
    BoxSize = header.BoxSize;
    RedShift = header.redshift;
    sprintf( LogBuf, "NumPart = %ld, N_Gas = %ld", NumPart, N_Gas );
    print_log( LogBuf );
    allocate_memory();
    show_header( header );
    for ( blk=0; blk<IO_NBLOCKS; blk++ ) {
        for ( pt=0, offset=0; pt<6; pt++ ) {
            for (file=0; file<All.NumFiles; file++) {
                if ( All.NumFiles < 2 )
                    sprintf( file_name, "%s_%03d.hdf5", All.FilePrefix, All.StartSnapIndex + ThisTask );
                else
                    sprintf( file_name, "%s_%03d.%3i.hdf5",
                            All.FilePrefix, All.StartSnapIndex + ThisTask, file );
                    read_header( file_name );
                    if ( blockpresent( blk, pt ) ) {
                        nbytes = get_block_nbytes( blk );
                        if ( All.BufferSize * 1024 * 1024 < nbytes * header.npart[pt] ){
                            printf( "BufferSize is too small.\n" );
                            endrun( 1 );
                        }
                        hdf5_file = H5Fopen( file_name, H5F_ACC_RDWR, H5P_DEFAULT );
                        get_dataset_name( blk, buf );
                        get_hdf5_native_type( blk, &hdf5_type );
                        sprintf( LogBuf, "[%i] reading %s ...", pt, buf );
                        print_log( LogBuf );
                        sprintf( buf1, "PartType%i/%s", pt, buf );
                        hdf5_dataset = H5Dopen( hdf5_file, buf1 );
                        herr = H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, CommBuffer );
                        empty_buffer( blk, offset, pt );
                        offset += header.npart[pt];
                        H5Dclose( hdf5_dataset );
                        H5Tclose( hdf5_type );
                        H5Fclose( hdf5_file );
                }
            }
        }
    }
    for ( pt=0, offset=0; pt<6; pt++ ) {
        num = header.npartTotal[pt] + ( ( ( long long )header.npartTotalHighWord[pt] ) << 32 );
        if ( num != 0 && header.mass[pt] != 0 ){
            sprintf( LogBuf, "copy particle %d mass from header to paritcle data", pt );
            print_log( LogBuf );
            for ( i=offset; i<offset+num; i++ )
                P[i].Mass = header.mass[pt];
        }
        offset += num;
    }
    print_log( "read data ... done. " );
    print_log( sep_str );
    find_id();
}

