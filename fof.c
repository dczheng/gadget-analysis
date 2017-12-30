#include "allvars.h"

long npart, offset;
double LinkL, rhodm;

void fof_allocate() {
    size_t bytes;

    if ( !(fof_info = malloc( bytes = npart * sizeof( struct fof_info_struct ) ) ) )
        printf( "failed to allocate memory for `fof_info` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `fof_info` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(fof_Next = malloc( bytes = npart * sizeof( long ) ) ) )
        printf( "failed to allocate memory for `fof_Next` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `fof_Next` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(Ngblist = malloc( bytes = npart * sizeof( long ) ) ) )
        printf( "failed to allocate memory for `Ngblist` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `Ngblist` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    print_log( sep_str );
}

void fof_free() {
    free( fof_info );
    free( fof_Next );
    free( Ngblist );
}

int fof_compare_len( const void *a, const void *b ) {
    return ( ( ( ( struct fof_info_struct *)a )->Len < ( ( struct fof_info_struct * )b )->Len ) ? 1 : -1 );
}

void fof_find_groups( int pt ) {
    int ngbnum, k, ngbmax;
    double pos[3];
    long i, ii, p, s, j, ss;
    //int flag=0;
    ngbmax = 0;
    print_log( "start fof find groups ..." );
    for ( i=offset; i<offset+npart; i++ ) {
        ii = i-offset;
        if ( !P[i].Flag ) continue;
        for( j=0; j<3; j++ )
            pos[j] = P[i].Pos[j];
        ngbnum = ngb_fof( pos, LinkL, pt );
        //printf( "%i\n", ngbnum );
        //if ( flag++ > 30 ) endrun( 20171230 );
        if ( ngbnum > ngbmax ) ngbmax = ngbnum;
        for ( k=0; k<ngbnum; k++ ) {
            j = Ngblist[k];
            if ( fof_info[ii].Head != fof_info[j].Head ) {
                if ( fof_info[fof_info[ii].Head].Len > fof_info[fof_info[j].Head].Len ) {
                    p = ii;
                    s = j;
                }
                else{
                    p = j;
                    s = ii;
                }
                fof_Next[ fof_info[fof_info[p].Head].Tail ] = fof_info[s].Head;
                fof_info[ fof_info[p].Head ].Tail = fof_info[ fof_info[s].Head ].Tail;
                fof_info[ fof_info[p].Head ].Len += fof_info[ fof_info[s].Head ].Len;
                fof_info[ fof_info[s].Head ].Len = 0;
                ss = fof_info[s].Head;
                do
                    fof_info[ss].Head = fof_info[p].Head;
                while( (ss=fof_Next[ss]) >= 0 );
            }
        }
    }
    qsort( fof_info, npart, sizeof( struct fof_info_struct ), fof_compare_len );
    sprintf( LogBuf, "the maximum number of ngb: %i", ngbmax );
    print_log( LogBuf );
    print_log( "fof find groups ... done" );
}

void fof_compute_group_properties( int pt ) {
    int i, num, j, k;
    long ii, p;
    double mass;
    print_log( "fof compute groups properties ... " );
    Ngroups = 0;
    for ( i=0; i<npart; i++ ) {
        if ( fof_info[i].Len < All.FofMinLen ) break;
        p = fof_info[i].Head;
        fof_info[i].mass = 0;
        for ( k=0; k<3; k++ ){
            fof_info[i].cm[k] = 0;
            fof_info[i].vel[k] = 0;
        }
        fof_info[i].vr200 = 0;
        for ( j=0; j<fof_info[i].Len; j++ ) {
            ii = p + offset;
            p = fof_Next[p];
            mass = ( header.mass[pt] != 0 ) ? header.mass[pt] : P[ii].Mass;
            fof_info[i].mass += mass;
            for ( k=0; k<3; k++ ){
                fof_info[i].cm[k] += mass * P[ii].Pos[k];
                fof_info[i].vel[k] += mass * P[ii].Vel[k];
            }
        }
        for ( k=0; k<3; k++ ){
            fof_info[i].cm[k] /= fof_info[i].mass;
            fof_info[i].vel[k] /= fof_info[i].mass;
        }
        fof_info[i].vr200 = pow( fof_info[i].mass /
            ( All.CriticalDensity * 200 * 4.0 / 3.0 * PI ), 1.0/3.0 );
    }
    Ngroups = i;
    sprintf( LogBuf, "The number of groups with at least %i particles: %li", All.FofMinLen, Ngroups );
    print_log( LogBuf );
    sprintf( LogBuf, "Largest group has %li particles", fof_info[0].Len );
    print_log( LogBuf );
    print_log( "fof compute groups properties ... done " );
}

void fof_test() {
    int i, index;
    long p, j;
    FILE *fd;
    index = 0;
    fd = fopen( "fof.txt", "w" );
    p = fof_info[index].Head;
    for ( i=0; i<fof_info[index].Len; i++ ){
        j = p + offset;
        fprintf( fd, "%g %g %g\n",
                P[j].Pos[0],
                P[j].Pos[1],
                P[j].Pos[2] );
        p = fof_Next[p];
    }
    fclose( fd );
}

void fof_save_groups() {
    hid_t hdf5_file, hdf5_dataset, hdf5_dataspace, hdf5_attribute, hdf5_type;
    herr_t herr;
    long *buf1, i, j;
    double *buf2;
    int ndims;
    hsize_t dims[2];
    print_log( "fof save groups ... " );
    sprintf( All.FofFileName, "%s.hdf5", All.FofFileName );
    hdf5_file = H5Fcreate( All.FofFileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_file, "GroupsNumberAboveMinLength", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &Ngroups );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_file, "MinLength", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &All.FofMinLen );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    ndims = 1;
    dims[0] = npart;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    hdf5_type = H5Tcopy( H5T_NATIVE_UINT64 );
    hdf5_dataset = H5Dcreate( hdf5_file, "Next", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, fof_Next );
    H5Dclose( hdf5_dataset );
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );
    /*************************************/
    if( !(buf1 = malloc( sizeof( long ) * Ngroups )) ){
        printf( "failed allocate `buf1` to save groups !\n" );
        fof_info[i].Head;
        endrun( 20171229 );
    }
    ndims = 1;
    dims[0] = Ngroups;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    hdf5_type = H5Tcopy( H5T_NATIVE_UINT64 );
    /*************************************/
    for ( i=0; i<Ngroups; i++ ) {
        buf1[i] = fof_info[i].Head;
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "Head", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf1 );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ ) {
        buf1[i] = fof_info[i].Len;
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "Length", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf1 );
    H5Dclose( hdf5_dataset );
    /*************************************/
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );
    free( buf1 );
    /*************************************/
    ndims = 1;
    dims[0] = Ngroups;
    if( !(buf2 = malloc( sizeof( double ) * Ngroups * 3 )) ){
        printf( "failed allocate `buf2` to save groups !\n" );
        fof_info[i].Head;
        endrun( 20171229 );
    }
    hdf5_type = H5Tcopy( H5T_NATIVE_DOUBLE );
    /*************************************/
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    for ( i=0; i<Ngroups; i++ ) {
        buf2[i] = fof_info[i].mass;
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "Mass", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf2 );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ ) {
        buf2[i] = fof_info[i].vr200;
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "VirialR200", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf2 );
    H5Dclose( hdf5_dataset );
    H5Sclose( hdf5_dataspace );
    /*************************************/
    ndims = 2;
    dims[0] = Ngroups;
    dims[1] = 3;
    /*************************************/
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    for ( i=0; i<Ngroups; i++ ) {
        for ( j=0; j<3; j++ )
        buf2[i*3+j] = fof_info[i].cm[j];
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "CenterOfMass", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf2 );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ ) {
        for ( j=0; j<3; j++ )
            buf2[i*3+j] = fof_info[i].vel[j];
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "VelocityOfCenterOfMass", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf2 );
    H5Dclose( hdf5_dataset );
    H5Sclose( hdf5_dataspace );
    /*************************************/
    H5Tclose( hdf5_type );
    free( buf2 );
    /*************************************/

    H5Fclose( hdf5_file );
    print_log( "fof save groups ... done" );
}

void fof( int pt ) {
    long i;
    double masstot, mass;
    print_log( "start fof ..." );
    npart = get_particle_num( pt );
    offset = get_particle_offset( pt );
    fof_allocate();
    tree_build( pt );
    for ( i=0; i<npart; i++ ) {
        fof_info[i].Head = fof_info[i].Tail = i;
        fof_info[i].Len = 1;
        fof_Next[i] = -1;
    }
    rhodm = (All.Omega0-All.OmegaBaryon) * 3 * SQR( All.Hubble ) / ( 8 * PI * All.G );
    if ( header.mass[pt] != 0 )
        mass = header.mass[pt];
    else{
        for ( i=offset, masstot=0; i<offset+npart; i++ )
            masstot += P[i].Mass;
        mass = masstot / npart;
    }
    LinkL = All.LinkLength * pow( mass / rhodm, 1.0/3 );
    sprintf( LogBuf, "critical density of dark matter: %g\n"
            "comoving linking lenght %g", rhodm, LinkL );
    print_log( LogBuf );
    fof_find_groups( pt );
    fof_compute_group_properties( pt );
    //fof_test();
    tree_free();
    print_log( "fof ... done" );
    print_log( sep_str );
}
