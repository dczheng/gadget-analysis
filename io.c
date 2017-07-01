#include "allvars.h"

void test();

void read_header( char *fn ) {

    //fprintf( stdout, "reading header From %s\n", file_name );
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

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_DoublePrecision");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision);
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
        fprintf( stdout, "%li ", header.npart[i] );
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
    fprintf( stdout, "%-25s: %i\n", "flag_sfr",               header.flag_sfr );
    fprintf( stdout, "%-25s: %i\n", "flag_feedback",          header.flag_feedback );
    fprintf( stdout, "%-25s: %i\n", "num_files",              header.num_files );
    fprintf( stdout, "%-25s: %lf\n", "BoxSize",               header.BoxSize );
    fprintf( stdout, "%-25s: %lf\n", "Omega0",                header.Omega0 );
    fprintf( stdout, "%-25s: %lf\n", "OmegaLambda",           header.OmegaLambda );
    fprintf( stdout, "%-25s: %i\n", "flag_stellarge",         header.flag_stellarage );
    fprintf( stdout, "%-25s: %i\n", "flag_metals",            header.flag_metals );
    fprintf( stdout, "%-25s: %i\n", "flag_entropy_instead_u", header.flag_entropy_instead_u );
    fprintf( stdout, "%-25s: %i\n", "flag_doubleprecision",   header.flag_doubleprecision );
    fprintf( stdout, "%-25s: %i\n", "flag_ic_info",           header.flag_ic_info );
    fprintf( stdout, "%-25s: %f\n", "lpt_scalingfactor",      header.lpt_scalingfactor );
    fputs( sep_str, stdout );
}

void write_header( char *fn, struct io_header header ) {

    fprintf(  stdout, "write header To %s\n", fn  );
    hdf5_file = H5Fcreate( fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    //hdf5_file = H5Fopen( fn, H5F_ACC_RDWR, H5P_DEFAULT );
    hdf5_group = H5Gcreate( hdf5_file, "/header", 0 );
    hsize_t adim[1] = { 6 };

    hdf5_dataspace = H5Screate( H5S_SIMPLE );
    H5Sset_extent_simple( hdf5_dataspace, 1, adim, NULL );
    hdf5_attribute = H5Acreate( hdf5_group, "NumPart_Thisfile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
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
    hdf5_attribute = H5Acreate( hdf5_group, "NumfilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
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

void write_file( char *fn, struct io_header header, struct Particle_Struct *Particle ) {
    write_header( fn, header );
    long pt;
    char buf[50];
    hdf5_file = H5Fopen( fn, H5F_ACC_RDWR, H5P_DEFAULT );
    for ( pt=0; pt<2; pt++ ) {

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
                    H5P_DEFAULT, Particle[pt].pos );
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
        hdf5_type = H5Tcopy( H5T_NATIVE_FLOAT );
        hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
        hdf5_dataset = H5Dcreate( hdf5_group, "Masses", hdf5_type,
                    hdf5_dataspace, H5P_DEFAULT );
        H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL,
                    H5P_DEFAULT, Particle[pt].m );
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

        if ( 0 == pt ) {
            ndims = 1;
            dims[0] = header.npart[pt];
            hdf5_type = H5Tcopy( H5T_NATIVE_FLOAT );
            hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
            hdf5_dataset = H5Dcreate( hdf5_group, "InternalEnergy", hdf5_type,
                    hdf5_dataspace, H5P_DEFAULT );
            H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL,
                    H5P_DEFAULT, Particle[pt].u );
            H5Tclose( hdf5_type );
            H5Sclose( hdf5_dataspace );
        }

        H5Gclose( hdf5_group );
    }
    H5Fclose( hdf5_file );
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
        case IO_MASS:
        case IO_U:
        case IO_RHO:
        case IO_POT:
        case IO_ELEC:
            block_nbytes = sizeof( float );
            break;
        case IO_ID:
            block_nbytes = sizeof( MyIDType );
            break;
    }
    return block_nbytes;
}

void get_block_dims( int pt, enum iofields blk, hsize_t *dims ) {
    int dim2 = 0;
    switch ( blk ) {
        case IO_POS:
        case IO_VEL:
        case IO_ACCEL:
        case IO_MAG:
            dims[0] = header.npart[pt];
            dims[1] = 3;
            break;
        case IO_MASS:
        case IO_U:
        case IO_RHO:
        case IO_POT:
        case IO_ELEC:
        case IO_ID:
            dims[0] = header.npart[pt];
            dims[1] = 1;
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
        case IO_ELEC:
            strcpy( buf, "ElectronAbundance" );
            break;
        case IO_ID:
            strcpy( buf, "ParticleIDs" );
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
            *hdf5_type = H5Tcopy(H5T_NATIVE_FLOAT);
            break;
    }
}

void read_block( int pt, void *data, enum iofields blk ) {
    char buf[50], buf1[20], file_name[FILENAME_MAX];
    void *p;
    long i;
    int nbytes;
    get_hdf5_native_type( blk, &hdf5_type );
    get_dataset_name( blk, buf1 );
    nbytes = get_block_nbytes( blk );
    sprintf( buf, "PartType%i/%s", pt, buf1 );
    p = data;
    for ( i=0; i<header.num_files; i++ ) {
        sprintf( file_name, "%s.%i.hdf5", file_Prefix, i );
        if ( header.num_files < 2 )
            sprintf( file_name, "%s.hdf5", file_Prefix);
        read_header( file_name );
        hdf5_file = H5Fopen( file_name , H5F_ACC_RDWR, H5P_DEFAULT );
        hdf5_dataset = H5Dopen( hdf5_file, buf );
        herr = H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, p );
        p = p + header.npart[pt]*nbytes;
        H5Dclose( hdf5_dataset );
        H5Fclose( hdf5_file );
    }
    H5Tclose( hdf5_type );
}

void read_all_data() {
    int pt, blk, nbytes, dim2;
    long i;
    char file_name[FILENAME_MAX], buf[20];
    fputs( sep_str, stdout );
    fputs( "read all data ...\n", stdout );
    sprintf( file_name, "%s.%i.hdf5", file_Prefix, 0 );
    if ( Num_files < 2 )
        sprintf( file_name, "%s.hdf5", file_Prefix );
    read_header( file_name );
    show_header( header );
    for ( pt=0; pt<6; pt++ ) {
        Particle[pt].num = header.npartTotalHighWord[pt];
        Particle[pt].num = header.npartTotal[pt] + ( Particle[pt].num<<32 );
        //fprintf( stdout, "%li\n", Particle[pt].num );
        fputs( sep_str, stdout );
        if ( Particle[pt].num != 0 ){
            fprintf( stdout, "Particle %i\n", pt );
            for ( blk=0; blk<IO_NBLOCKS; blk++ ){
                switch ( blk ){
                    case IO_POS:
                        nbytes = get_block_nbytes( blk );
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "reading %s\n", buf);
                        Particle[pt].pos = ( float* ) malloc( nbytes * Particle[pt].num );
                        read_block( pt, (void*)(Particle[pt].pos), blk );
                        break;
                    case IO_VEL:
                        nbytes = get_block_nbytes( blk );
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "reading %s\n", buf);
                        Particle[pt].vel = ( float* ) malloc( nbytes * Particle[pt].num );
                        read_block( pt, (void*)(Particle[pt].vel), blk );
                        break;
                    case IO_ACCEL:
                        nbytes = get_block_nbytes( blk );
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "reading %s\n", buf);
                        Particle[pt].accel = ( float* ) malloc( nbytes * Particle[pt].num );
                        read_block( pt, (void*)(Particle[pt].accel), blk );
                        break;
                    case IO_MAG:
                        if ( pt>0 ) break;
                        nbytes = get_block_nbytes( blk );
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "reading %s\n", buf);
                        Particle[pt].mag = ( float* ) malloc(  nbytes * Particle[pt].num );
                        read_block( pt, (void*)(Particle[pt].mag), blk );
                        break;
                    case IO_MASS:
                        nbytes = get_block_nbytes( blk );
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "reading %s\n", buf);
                        Particle[pt].m = ( float* ) malloc( nbytes * Particle[pt].num );
                        if ( header.mass[pt] == 0 )
                            read_block( pt, (void*)(Particle[pt].m), blk );
                        else
                            for ( i=0; i<Particle[pt].num; i++ )
                                Particle[pt].m[i] = header.mass[pt];
                        break;
                    case IO_U:
                        if ( pt>0 ) break;
                        nbytes = get_block_nbytes( blk );
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "reading %s\n", buf);
                        Particle[pt].u = ( float* ) malloc( nbytes * Particle[pt].num );
                        read_block( pt, (void*)(Particle[pt].u), blk );
                        break;
                    case IO_POT:
                        nbytes = get_block_nbytes( blk );
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "reading %s\n", buf);
                        Particle[pt].pot = ( float* ) malloc( nbytes * Particle[pt].num );
                        read_block( pt, (void*)(Particle[pt].pot), blk );
                        break;
                    case IO_ELEC:
                        if ( pt>0 ) break;
                        nbytes = get_block_nbytes( blk );
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "reading %s\n", buf);
                        Particle[pt].elec = ( float* ) malloc( nbytes * Particle[pt].num );
                        read_block( pt, (void*)(Particle[pt].elec), blk );
                        break;
                    case IO_RHO:
                        if ( pt>0 ) break;
                        nbytes = get_block_nbytes( blk );
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "reading %s\n", buf);
                        Particle[pt].rho = ( float* ) malloc( nbytes * Particle[pt].num );
                        read_block( pt, (void*)(Particle[pt].rho), blk );
                        break;
                    case IO_ID:
                        nbytes = get_block_nbytes( blk );
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "reading %s\n", buf);
                        Particle[pt].id = ( MyIDType* ) malloc( nbytes * Particle[pt].num );
                        read_block( pt, (void*)(Particle[pt].id), blk );
                        break;
                }
            }
        }
        fputs( sep_str, stdout );
    }
    fputs( sep_str, stdout );
    //test();
}

void free_all_memory() {
    int pt, blk;
    char buf[20];
    fputs( sep_str, stdout );
    fputs( "free all memory ...\n", stdout );
    for ( pt=0; pt<6; pt++) {
        fputs( sep_str, stdout );
        if ( Particle[pt].num != 0 ){
            fprintf( stdout, "Particle %i\n", pt );
            for ( blk=0; blk<IO_NBLOCKS; blk++ ){
                switch ( blk ){
                    case IO_POS:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "free memory: %s\n", buf);
                        free( Particle[pt].pos );
                        break;
                    case IO_VEL:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "free memory: %s\n", buf);
                        free( Particle[pt].vel );
                        break;
                    case IO_ACCEL:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "free memory: %s\n", buf);
                        free( Particle[pt].accel );
                        break;
                    case IO_MAG:
                        if ( pt>0 ) break;
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "free memory: %s\n", buf);
                        free( Particle[pt].mag );
                        break;
                    case IO_MASS:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "free memory: %s\n", buf);
                        free( Particle[pt].m );
                        break;
                    case IO_U:
                        if ( pt>0 ) break;
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "free memory: %s\n", buf);
                        free( Particle[pt].u );
                        break;
                    case IO_POT:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "free memory: %s\n", buf);
                        free( Particle[pt].pot );
                        break;
                    case IO_ELEC:
                        if ( pt>0 ) break;
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "free memory: %s\n", buf);
                        free( Particle[pt].elec );
                        break;
                    case IO_RHO:
                        if ( pt>0 ) break;
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "free memory: %s\n", buf);
                        free( Particle[pt].rho );
                        break;
                    case IO_ID:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "free memory: %s\n", buf);
                        free( Particle[pt].id );
                        break;
                }
            }
        }
        fputs( sep_str, stdout );
    }
    fputs( sep_str, stdout );
}

void test() {
    int pt, blk, i;
    char buf[20];
    int test_num;
    fputs( sep_str, stdout );
    fputs( "test ...\n", stdout );
    test_num = 3;
    for ( pt=0; pt<6; pt++) {
        if ( pt>=2 ) continue;
        fputs( sep_str, stdout );
        fprintf( stdout, "Particle %i\n", pt );
        if ( Particle[pt].num != 0 )
            for ( blk=0; blk<IO_NBLOCKS; blk++ ){
                switch ( blk ){
                    case IO_POS:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "%s\n", buf);
                        for ( i=0; i<test_num; i++ ) {
                            fprintf( stdout, "%10.4f %10.4f %10.4f\n",
                                    Particle[pt].pos[i*3+0],
                                    Particle[pt].pos[i*3+1],
                                    Particle[pt].pos[i*3+2] );
                        }
                        break;
                    case IO_VEL:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "%s\n", buf);
                        for ( i=0; i<test_num; i++ ) {
                            fprintf( stdout, "%10.4f %10.4f %10.4f\n",
                                    Particle[pt].vel[i*3+0],
                                    Particle[pt].vel[i*3+1],
                                    Particle[pt].vel[i*3+2] );
                        }
                        break;
                    case IO_ACCEL:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "%s\n", buf);
                        for ( i=0; i<test_num; i++ ) {
                            fprintf( stdout, "%10.4f %10.4f %10.4f\n",
                                    Particle[pt].accel[i*3+0],
                                    Particle[pt].accel[i*3+1],
                                    Particle[pt].accel[i*3+2] );
                        }
                        break;
                    case IO_MAG:
                        if ( pt>0 ) break;
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "%s\n", buf);
                        for ( i=0; i<test_num; i++ ) {
                            fprintf( stdout, "%e %e %e\n",
                                    Particle[pt].mag[i*3+0],
                                    Particle[pt].mag[i*3+1],
                                    Particle[pt].mag[i*3+2] );
                        }
                        break;
                    case IO_MASS:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "%s\n", buf);
                        for ( i=0; i<test_num; i++ ) {
                            fprintf( stdout, "%10.4f\n", Particle[pt].m[i] );
                        }
                        break;
                    case IO_U:
                        if ( pt>0 ) break;
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "%s\n", buf);
                        for ( i=0; i<test_num; i++ ) {
                            fprintf( stdout, "%e\n", Particle[pt].u[i] );
                        }
                        break;
                    case IO_POT:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "%s\n", buf);
                        for ( i=0; i<test_num; i++ ) {
                            fprintf( stdout, "%e\n", Particle[pt].pot[i] );
                        }
                        break;
                    case IO_ELEC:
                        if ( pt>0 ) break;
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "%s\n", buf);
                        for ( i=0; i<test_num; i++ ) {
                            fprintf( stdout, "%e\n", Particle[pt].elec[i] );
                        }
                        break;
                    case IO_RHO:
                        if ( pt>0 ) break;
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "%s\n", buf);
                        for ( i=0; i<test_num; i++ ) {
                            fprintf( stdout, "%e\n", Particle[pt].rho[i] );
                        }
                        break;
                    case IO_ID:
                        get_dataset_name( blk, buf );
                        fprintf( stdout, "%s\n", buf);
                        for ( i=0; i<test_num; i++ ) {
#ifdef LONGIDS
                            fprintf( stdout, "%llu\n", Particle[pt].id[i] );
#else
                            fprintf( stdout, "%lu\n", Particle[pt].id[i] );
#endif
                        }
                        break;
                }
            }
        fputs( sep_str, stdout );
    }
    fputs( sep_str, stdout );
}

