#include "allvars.h"

double LinkL, rhodm;

void fof_allocate() {
    size_t bytes;

    mymalloc( FoFProps, NumPart * sizeof( struct fof_properties ) );
    mymalloc( FoFNext, NumPart * sizeof( long ) );
    mymalloc( Ngblist, NumPart * sizeof( long ) );

    writelog( sep_str );
}

void fof_free() {
    myfree( FoFProps );
    myfree( FoFNext );
    myfree( Ngblist );
}

int fof_compare_len( const void *a, const void *b ) {
    return ( ( ( ( struct fof_properties *)a )->Len < ( ( struct fof_properties * )b )->Len ) ? 1 : -1 );
}

void fof_find_groups() {
    int ngbnum, k, ngbmax;
    double pos[3];
    long i, p, s, j, ss,
        *Head, *Tail,*Len;
    time_start();
    ngbmax = 0;
    writelog( sep_str );
    writelog( "Start FoF find groups ...\n" );

    mymalloc( Head, NumPart * sizeof( long ) );
    mymalloc( Tail, NumPart * sizeof( long ) );
    mymalloc( Len, NumPart * sizeof( long ) );

    for ( i=0; i<NumPart; i++ ) {
        Head[i] = Tail[i] = i;
        Len[i] = 1;
        FoFNext[i] = -1;
    }

    for ( i=0; i<NumPart; i++ ) {
        if ( !( ( 1 << P[i].Type ) & All.TreePartType ) )
            continue;
        for( j=0; j<3; j++ )
            pos[j] = P[i].Pos[j];
        ngbnum = ngb_fof( pos, LinkL );
        if ( ngbnum > ngbmax ) ngbmax = ngbnum;
        ngbmax = ( ngbnum > ngbmax ) ? ngbnum : ngbmax;
        for ( k=0; k<ngbnum; k++ ) {
            j = Ngblist[k];
            if ( Head[i] != Head[j] ) {
                if ( Len[Head[i]] > Len[Head[j]] ) {
                    p = i;
                    s = j;
                }
                else{
                    p = j;
                    s = i;
                }
                FoFNext[Tail[Head[p]]] = Head[s];
                Tail[Head[p]] = Tail[Head[s]];
                Len[Head[p]] += Len[Head[s]];
                Len[Head[s]] = 1;

                ss = Head[s];
                do
                    Head[ss] = Head[p];
                while( (ss=FoFNext[ss]) >= 0 );
            }
        }
    }

    for ( i=0; i<NumPart; i++ ) {
        FoFProps[i].Head = Head[i];
        FoFProps[i].Tail = Tail[i];
        FoFProps[i].Len = Len[i];
    }

    myfree( Head );
    myfree( Tail );
    myfree( Len );

    writelog( "the maximum number of ngb: %i\n", ngbmax );
    writelog( "FoF find groups ... done\n" );
    time_end();
    writelog( sep_str );
}

void fof_compute_group_properties() {
    long p, i, j, k, p0;

    time_start();
    writelog( "FoF compute groups properties ... \n" );

    Ngroups = 0;
    for ( i=0; i<NumPart; i++ ) {
        if ( FoFProps[i].Len >= All.FofMinLen ) {
            FoFProps[ Ngroups++ ] = FoFProps[i];
        }
    }

    qsort( FoFProps, Ngroups, sizeof( struct fof_properties ), fof_compare_len );

    writelog( "The number of groups with at least %i particles: %li\n", All.FofMinLen, Ngroups );
    writelog( "Largest group has %li particles\n", FoFProps[0].Len );

    for ( i=0; i<Ngroups; i++ ) {
        p0 = p = FoFProps[i].Head;
        FoFProps[i].mass = 0;
        for ( k=0; k<3; k++ ){
            FoFProps[i].cm[k] = 0;
            FoFProps[i].vel[k] = 0;
        }
        FoFProps[i].vr200 = 0;
        for ( j=0; j<FoFProps[i].Len; j++ ) {
            FoFProps[i].mass += P[p].Mass;
            for ( k=0; k<3; k++ ){
                FoFProps[i].cm[k] += P[p].Mass *
                    ( PERIODIC( P[p].Pos[k]-P[p0].Pos[k] ) + P[p0].Pos[k] );
                FoFProps[i].vel[k] += P[p].Mass * P[p].Vel[k];
            }
            p = FoFNext[p];
        }
        if ( FoFProps[i].mass == 0 ) {
            printf( "FoFProps[%li] is zeros !!!\n", i );
            endrun( 20180507 );
        }
        for ( k=0; k<3; k++ ){
            FoFProps[i].cm[k] /= FoFProps[i].mass;
            FoFProps[i].vel[k] /= FoFProps[i].mass;
        }
        FoFProps[i].vr200 = pow( FoFProps[i].mass /
            ( All.RhoCrit * 200 * 4.0 / 3.0 * PI ), 1.0/3.0 );
    }
    writelog( "FoF compute groups properties ... done\n" );
    time_end();
    writelog( sep_str );

}

void fof_test() {
    int i, index;
    long p;
    FILE *fd;
    index = 0;
    fd = fopen( "fof.txt", "w" );
    p = FoFProps[index].Head;
    for ( i=0; i<FoFProps[index].Len; i++ ){
        fprintf( fd, "%g %g %g\n",
                P[p].Pos[0],
                P[p].Pos[1],
                P[p].Pos[2] );
        p = FoFNext[p];
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
    writelog( "FoF save groups ...\n" );
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
    dims[0] = NumPart;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    hdf5_type = H5Tcopy( H5T_NATIVE_UINT64 );
    hdf5_dataset = H5Dcreate( hdf5_file, "Next", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, FoFNext );
    H5Dclose( hdf5_dataset );
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );
    /*************************************/

    mymalloc( buf1, Ngroups * sizeof( long ) );

    ndims = 1;
    dims[0] = Ngroups;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    hdf5_type = H5Tcopy( H5T_NATIVE_UINT64 );
    /*************************************/
    for ( i=0; i<Ngroups; i++ ) {
        buf1[i] = FoFProps[i].Head;
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "Head", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf1 );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ ) {
        buf1[i] = FoFProps[i].Len;
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "Length", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf1 );
    H5Dclose( hdf5_dataset );
    /*************************************/
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );
    myfree( buf1 );
    /*************************************/
    ndims = 1;
    dims[0] = Ngroups;

    mymalloc( buf2, Ngroups * sizeof( double ) * 3 );

    hdf5_type = H5Tcopy( H5T_NATIVE_DOUBLE );
    /*************************************/
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    for ( i=0; i<Ngroups; i++ ) {
        buf2[i] = FoFProps[i].mass;
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "Mass", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf2 );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ ) {
        buf2[i] = FoFProps[i].vr200;
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
        buf2[i*3+j] = FoFProps[i].cm[j];
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "CenterOfMass", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf2 );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroups; i++ ) {
        for ( j=0; j<3; j++ )
            buf2[i*3+j] = FoFProps[i].vel[j];
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "VelocityOfCenterOfMass", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf2 );
    H5Dclose( hdf5_dataset );
    H5Sclose( hdf5_dataspace );
    /*************************************/
    H5Tclose( hdf5_type );
    myfree( buf2 );
    /*************************************/

    H5Fclose( hdf5_file );
    writelog( "FoF save groups ... done\n" );
}

void fof() {
    long i, npart;
    double masstot, mass;
    time_start();
    writelog( "Start FoF ...\n" );
    for ( i=0, npart=0, masstot=0; i<NumPart; i++ )
        if ( ( 1 << P[i].Type ) & All.TreePartType ) {
            npart ++;
            masstot += P[i].Mass;
        }
    writelog( "particle type: ` " );
    for ( i=0; i<6; i++ )
        if ( ( 1 << i ) & All.TreePartType )
            writelog( "%li ", i );
    writelog( "` used in FoF\n" )
    writelog( "total particle number is: %li\n", npart );

    fof_allocate();
    tree_build();

    rhodm = (All.Omega0-All.OmegaBaryon) * 3 *  SQR(All.Hubble) / ( 8 * PI * All.G );
    mass = masstot / npart;
    LinkL = All.LinkLength * pow( mass / rhodm, 1.0/3 );
    //LinkL = 65.5483;
    /*
    printf( "%.10f %.10f %.10f %.10f %.20f %.10f %.10f %li \n",
             All.Omega0, All.OmegaBaryon,
             All.Hubble, All.G, rhodm, M_PI, masstot, npart );
             */

    writelog( "critical density of dark matter: %g\n"
            "comoving linking lenght %g\n", rhodm, LinkL );
    fof_find_groups();
    fof_compute_group_properties();
    //fof_test();
    tree_free();
    writelog( "FoF ... done\n" );
    time_end();
    writelog( sep_str );
}
