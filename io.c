#include "allvars.h"

#if defined(RAD) 
void save_particle_radio() {

    char fn[ FILENAME_MAX ];
    int ndims;
    hsize_t dims[2];

    put_header( "save radio" );
    hid_t hdf5_file, hdf5_dataset, hdf5_dataspace, hdf5_attribute, hdf5_type;

    sprintf( fn, "%s/rad_%03i.hdf5", All.RadDir, SnapIndex );

    hdf5_file = H5Fcreate( fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_file, "NuNum", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &All.NuNum );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_file, "NuMin", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &All.NuMin );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_file, "NuMax", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &All.NuMax );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    ndims = 2;
    dims[0] = N_Gas;
    dims[1] = All.NuNum;

    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    hdf5_type = H5Tcopy( H5T_NATIVE_DOUBLE );
    hdf5_dataset = H5Dcreate( hdf5_file, "Radio", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, PartRad );
    H5Dclose( hdf5_dataset );
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );

    H5Fclose( hdf5_file );

}

int read_particle_radio() {

    double NuMin, NuMax;
    int nuN;
    char fn[ FILENAME_MAX ];
    hid_t hdf5_file, hdf5_dataset, hdf5_attribute, hdf5_type;

    writelog( "read radio ...\n" );

    sprintf( fn, "%s/rad_%03i.hdf5", All.RadDir, SnapIndex );

    hdf5_file = H5Fopen( fn, H5F_ACC_RDWR, H5P_DEFAULT );

    hdf5_attribute = H5Aopen_name( hdf5_file, "NuNum" );
    H5Aread( hdf5_attribute, H5T_NATIVE_INT, &nuN );
    H5Aclose( hdf5_attribute );

    if ( nuN != All.NuNum ) {
        H5Fclose( hdf5_file );
        return 0;
    }

    hdf5_attribute = H5Aopen_name( hdf5_file, "NuMin" );
    H5Aread( hdf5_attribute, H5T_NATIVE_DOUBLE, &NuMin );
    H5Aclose( hdf5_attribute );

    if ( NuMin != All.NuMin ) {
        H5Fclose( hdf5_file );
        return 0;
    }

    hdf5_attribute = H5Aopen_name( hdf5_file, "NuMax" );
    H5Aread( hdf5_attribute, H5T_NATIVE_DOUBLE, &NuMax );
    H5Aclose( hdf5_attribute );

    if ( NuMax != All.NuMax ) {
        H5Fclose( hdf5_file );
        return 0;
    }

    hdf5_type = H5Tcopy( H5T_NATIVE_DOUBLE );
    hdf5_dataset = H5Dopen( hdf5_file, "Radio" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, PartRad );
    H5Dclose( hdf5_dataset );
    H5Tclose( hdf5_type );

    H5Fclose( hdf5_file );

    return 1;

}
#endif

#ifdef FOF

void fof_save() {

    hid_t hdf5_file, hdf5_dataset, hdf5_dataspace, hdf5_attribute, hdf5_type;
    long *buf1, i, j;
    double *buf2;
    char *buf;
    int ndims;
    char fn[50];
    hsize_t dims[2];
    writelog( "FoF save groups ...\n" );
    mytimer_start();

    sprintf( fn, "%s/fof_%03i.hdf5", All.FoFDir, SnapIndex );

    hdf5_file = H5Fcreate( fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_file, "GroupsNumberAboveMinLength", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &Ngroups );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_file, "MinLength", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &All.FoFMinLen );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    ndims = 1;
    dims[0] = NumPart;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    hdf5_type = H5Tcopy( H5T_NATIVE_INT64 );
    hdf5_dataset = H5Dcreate( hdf5_file, "Next", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, FoFNext );
    H5Dclose( hdf5_dataset );
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );


    if ( sizeof( double ) > sizeof( long ) ) {
        mymalloc1( buf, Ngroups * sizeof( double ) * 6 );
    }
    else {
        mymalloc1( buf, Ngroups * sizeof( long ) * 6 );
    }

    buf1 = (long *)buf;
    buf2 = (double *)buf;

    /*
    printf( "%i\n", Ngroups );
    endrun(20181107);
    */

    hdf5_type = H5Tcopy( H5T_NATIVE_UINT64 );
    /*****************int 1d********************/
    ndims = 1;
    dims[0] = Ngroups;

    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );

    for ( i=0; i<Ngroups; i++ ) {
        buf1[i] = Gprops[i].Head;
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "Head", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    for ( i=0; i<Ngroups; i++ ) {
        buf1[i] = Gprops[i].Len;
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "Length", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    H5Sclose( hdf5_dataspace );
    /*****************int 1d********************/

    /*****************int 2d********************/
    ndims = 2;
    dims[0] = Ngroups;
    dims[1] = 6;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );

    for ( i=0; i<Ngroups; i++ )
        for ( j=0; j<6; j++ )
            buf1[i*6+j] = Gprops[i].npart[j];
    hdf5_dataset = H5Dcreate( hdf5_file, "npart", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    H5Sclose( hdf5_dataspace );
    /*****************int 2d********************/
    H5Tclose( hdf5_type );

    hdf5_type = H5Tcopy( H5T_NATIVE_DOUBLE );
    /*****************real 1d********************/
    ndims = 1;
    dims[0] = Ngroups;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );

    for ( i=0; i<Ngroups; i++ ) {
        buf2[i] = Gprops[i].mass;
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "Mass", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    for ( i=0; i<Ngroups; i++ ) {
        buf2[i] = Gprops[i].vr200;
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "VirialR200", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    for ( i=0; i<Ngroups; i++ ) {
        buf2[i] = Gprops[i].ek;
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "KineticEnergy", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    
    for ( i=0; i<Ngroups; i++ ) {
        buf2[i] = Gprops[i].v_mean;
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "VMEAN", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    for ( i=0; i<Ngroups; i++ ) {
        buf2[i] = Gprops[i].v_disp;
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "VDISP", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    H5Sclose( hdf5_dataspace );
    /*****************real 1d********************/

    /*****************real 2d********************/
    ndims = 2;
    dims[0] = Ngroups;
    dims[1] = 3;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );

    for ( i=0; i<Ngroups; i++ ) {
        for ( j=0; j<3; j++ ) {
            buf2[i*3+j] = Gprops[i].size[j];
        }
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "Size", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );


    for ( i=0; i<Ngroups; i++ )
        for ( j=0; j<3; j++ )
            buf2[i*3+j] = Gprops[i].cm[j];

    hdf5_dataset = H5Dcreate( hdf5_file, "CenterOfMass", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    for ( i=0; i<Ngroups; i++ )
        for ( j=0; j<3; j++ )
            buf2[i*3+j] = Gprops[i].vel[j];

    hdf5_dataset = H5Dcreate( hdf5_file, "CenterOfVelocity", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    H5Sclose( hdf5_dataspace );

    dims[1] = 6;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );

    for ( i=0; i<Ngroups; i++ )
        for ( j=0; j<6; j++ )
        buf2[i*6+j] = Gprops[i].mass_table[j];

    hdf5_dataset = H5Dcreate( hdf5_file, "MassTable", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    H5Sclose( hdf5_dataspace );

    /*****************real 2d********************/
    H5Tclose( hdf5_type );
    myfree( buf );

    H5Fclose( hdf5_file );
    mytimer_end();
    writelog( "FoF save groups ... done\n" );
}

void fof_read() {

    hid_t hdf5_file, hdf5_dataset, hdf5_attribute, hdf5_type;
    long *buf1, i, j;
    double *buf2;
    char *buf;
    char fn[ FILENAME_MAX ];
    writelog( "read fof...\n" );
    mytimer_start();

    sprintf( fn, "%s/fof_%03i.hdf5", All.FoFDir, SnapIndex );

    hdf5_file = H5Fopen( fn, H5F_ACC_RDWR, H5P_DEFAULT );

    hdf5_attribute = H5Aopen_name( hdf5_file, "GroupsNumberAboveMinLength" );
    H5Aread( hdf5_attribute, H5T_NATIVE_INT, &Ngroups );
    H5Aclose( hdf5_attribute );

    hdf5_attribute = H5Aopen_name( hdf5_file, "MinLength" );
    H5Aread( hdf5_attribute, H5T_NATIVE_INT, &All.FoFMinLen );
    H5Aclose( hdf5_attribute );

    mymalloc2( Gprops, Ngroups * sizeof( struct group_properties ) );
    mymalloc2( FoFNext, NumPart * sizeof( long ) );

    writelog( "Ngroups: %i, MinLength: %i\n", Ngroups, All.FoFMinLen );

    if ( sizeof( double ) > sizeof( long ) ) {
        mymalloc1( buf, Ngroups * sizeof( double ) * 6 );
    }
    else {
        mymalloc1( buf, Ngroups * sizeof( long ) * 6 );
    }

    buf1 = (long *)buf;
    buf2 = (double *)buf;


    /*****************int********************/
    hdf5_type = H5Tcopy( H5T_NATIVE_INT64 );
    writelog( "read Next ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "Next" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, FoFNext );
    H5Dclose( hdf5_dataset );
    H5Tclose( hdf5_type );

    hdf5_type = H5Tcopy( H5T_NATIVE_UINT64 );
    writelog( "read Head ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "Head" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ )
        Gprops[i].Head = buf1[i];

    writelog( "read Length ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "Length" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ )
        Gprops[i].Len = buf1[i];

    writelog( "read npart ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "npart" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ )
        for ( j=0; j<6; j++ )
            Gprops[i].npart[j] = buf1[ i*6+j ];

    /*****************int********************/
    H5Tclose( hdf5_type );

    hdf5_type = H5Tcopy( H5T_NATIVE_DOUBLE );
    /*****************real********************/

    writelog( "read Mass ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "Mass" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ )
        Gprops[i].mass = buf2[i];

    writelog( "read VirialR200 ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "VirialR200" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++)
        Gprops[i].vr200 = buf2[i];

    writelog( "read KineticEnergy ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "KineticEnergy" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++)
        Gprops[i].ek = buf2[i];

    writelog( "read v_mean ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "VMEAN" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++)
        Gprops[i].v_mean = buf2[i];

    writelog( "read v_disp ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "VDISP" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++)
        Gprops[i].v_disp = buf2[i];

    writelog( "read Size ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "Size" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++)
        for (j=0; j<3; j++)
            Gprops[i].size[j] = buf2[i*3+j];

    writelog( "read CenterOfMass ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "CenterOfMass" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ )
        for ( j=0; j<3; j++ )
            Gprops[i].cm[j] = buf2[i*3+j];

    writelog( "read CenterOfVelocity ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "CenterOfVelocity" );
     H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ )
        for ( j=0; j<3; j++ )
            Gprops[i].vel[j] = buf2[i*3+j];

    writelog( "read MassTable ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "MassTable" );
     H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ )
        for ( j=0; j<6; j++ )
            Gprops[i].mass_table[j] = buf2[ i*6+j ];

    /*****************real********************/
    H5Tclose( hdf5_type );

    myfree( buf );

    H5Fclose( hdf5_file );
    mytimer_end();
    writelog( "read fof... done.\n" );
}
#endif

