#include "allvars.h"

#ifdef FOF
double LinkL, rhodm;

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
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &Ngroup );
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
        mymalloc1( buf, Ngroup * sizeof( double ) * 6 );
    }
    else {
        mymalloc1( buf, Ngroup * sizeof( long ) * 6 );
    }

    buf1 = (long *)buf;
    buf2 = (double *)buf;

    /*
    printf( "%i\n", Ngroup );
    endrun(20181107);
    */

    hdf5_type = H5Tcopy( H5T_NATIVE_UINT64 );
    /*****************int 1d********************/
    ndims = 1;
    dims[0] = Ngroup;

    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );

    for ( i=0; i<Ngroup; i++ ) {
        buf1[i] = Gprops[i].Head;
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "Head", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    for ( i=0; i<Ngroup; i++ ) {
        buf1[i] = Gprops[i].Len;
    }
    hdf5_dataset = H5Dcreate( hdf5_file, "Length", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    H5Sclose( hdf5_dataspace );
    /*****************int 1d********************/

    /*****************int 2d********************/
    ndims = 2;
    dims[0] = Ngroup;
    dims[1] = 6;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );

    for ( i=0; i<Ngroup; i++ )
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
    dims[0] = Ngroup;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );

    for ( i=0; i<Ngroup; i++ ) {
        buf2[i] = Gprops[i].mass;
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "Mass", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    for ( i=0; i<Ngroup; i++ ) {
        buf2[i] = Gprops[i].vr200;
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "VirialR200", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    for ( i=0; i<Ngroup; i++ ) {
        buf2[i] = Gprops[i].ek;
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "KineticEnergy", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    
    for ( i=0; i<Ngroup; i++ ) {
        buf2[i] = Gprops[i].v_mean;
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "VMEAN", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    for ( i=0; i<Ngroup; i++ ) {
        buf2[i] = Gprops[i].v_disp;
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "VDISP", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    H5Sclose( hdf5_dataspace );
    /*****************real 1d********************/

    /*****************real 2d********************/
    ndims = 2;
    dims[0] = Ngroup;
    dims[1] = 3;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );

    for ( i=0; i<Ngroup; i++ ) {
        for ( j=0; j<3; j++ ) {
            buf2[i*3+j] = Gprops[i].size[j];
        }
    }

    hdf5_dataset = H5Dcreate( hdf5_file, "Size", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );


    for ( i=0; i<Ngroup; i++ )
        for ( j=0; j<3; j++ )
            buf2[i*3+j] = Gprops[i].cm[j];

    hdf5_dataset = H5Dcreate( hdf5_file, "CenterOfMass", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    for ( i=0; i<Ngroup; i++ )
        for ( j=0; j<3; j++ )
            buf2[i*3+j] = Gprops[i].vel[j];

    hdf5_dataset = H5Dcreate( hdf5_file, "CenterOfVelocity", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );

    H5Sclose( hdf5_dataspace );

    dims[1] = 6;
    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );

    for ( i=0; i<Ngroup; i++ )
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
    H5Aread( hdf5_attribute, H5T_NATIVE_INT, &Ngroup );
    H5Aclose( hdf5_attribute );

    hdf5_attribute = H5Aopen_name( hdf5_file, "MinLength" );
    H5Aread( hdf5_attribute, H5T_NATIVE_INT, &All.FoFMinLen );
    H5Aclose( hdf5_attribute );

    mymalloc2( Gprops, Ngroup * sizeof( struct group_properties ) );
    mymalloc2( FoFNext, NumPart * sizeof( long ) );

    writelog( "Ngroup: %i, MinLength: %i\n", Ngroup, All.FoFMinLen );

    if ( sizeof( double ) > sizeof( long ) ) {
        mymalloc1( buf, Ngroup * sizeof( double ) * 6 );
    }
    else {
        mymalloc1( buf, Ngroup * sizeof( long ) * 6 );
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
    for ( i=0; i<Ngroup; i++ )
        Gprops[i].Head = buf1[i];

    writelog( "read Length ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "Length" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroup; i++ )
        Gprops[i].Len = buf1[i];

    writelog( "read npart ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "npart" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroup; i++ )
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
    for ( i=0; i<Ngroup; i++ )
        Gprops[i].mass = buf2[i];

    writelog( "read VirialR200 ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "VirialR200" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroup; i++)
        Gprops[i].vr200 = buf2[i];

    writelog( "read KineticEnergy ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "KineticEnergy" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroup; i++)
        Gprops[i].ek = buf2[i];

    writelog( "read v_mean ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "VMEAN" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroup; i++)
        Gprops[i].v_mean = buf2[i];

    writelog( "read v_disp ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "VDISP" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroup; i++)
        Gprops[i].v_disp = buf2[i];

    writelog( "read Size ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "Size" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroup; i++)
        for (j=0; j<3; j++)
            Gprops[i].size[j] = buf2[i*3+j];

    writelog( "read CenterOfMass ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "CenterOfMass" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroup; i++ )
        for ( j=0; j<3; j++ )
            Gprops[i].cm[j] = buf2[i*3+j];

    writelog( "read CenterOfVelocity ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "CenterOfVelocity" );
     H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroup; i++ )
        for ( j=0; j<3; j++ )
            Gprops[i].vel[j] = buf2[i*3+j];

    writelog( "read MassTable ...\n" );
    hdf5_dataset = H5Dopen( hdf5_file, "MassTable" );
     H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    for ( i=0; i<Ngroup; i++ )
        for ( j=0; j<6; j++ )
            Gprops[i].mass_table[j] = buf2[ i*6+j ];

    /*****************real********************/
    H5Tclose( hdf5_type );

    myfree( buf );

    H5Fclose( hdf5_file );
    mytimer_end();
    writelog( "read fof... done.\n" );
}

int fof_compare_len( const void *a, const void *b ) {
    return ( ( ( ( struct group_properties *)a )->Len < ( ( struct group_properties * )b )->Len ) ? 1 : -1 );
}

int fof_compare_mass( const void *a, const void *b ) {
    return ( ( ( ( struct group_properties *)a )->mass < ( ( struct group_properties * )b )->mass ) ? 1 : -1 );
}

void fof_find_groups() {

    int ngbnum, k, ngbmax;
    double pos[3];
    long i, p, s, j, ss,
        *Head, *Tail,*Len;
    ngbmax = 0;
    writelog( "Start FoF find groups ...\n" );
    mytimer_start();

    mymalloc1( Head, NumPart * sizeof( long ) );
    mymalloc1( Tail, NumPart * sizeof( long ) );
    mymalloc1( Len, NumPart * sizeof( long ) );
    mymalloc1( Ngblist, NumPart * sizeof( long ) );

    for ( i=0; i<NumPart; i++ ) {
        Head[i] = Tail[i] = i;
        Len[i] = 1;
        FoFNext[i] = -1;
    }

    for ( i=0; i<NumPart; i++ ) {

        //printf( "%li\n", i );

        if ( !( ( 1 << P[i].Type ) & All.TreePartType ) )
            continue;

        for( j=0; j<3; j++ )
            pos[j] = P[i].Pos[j];

        ngbnum = ngb_fof( pos, LinkL );
        if ( ngbnum > ngbmax ) ngbmax = ngbnum;

        ngbmax = ( ngbnum > ngbmax ) ? ngbnum : ngbmax;

        for ( k=0; k<ngbnum; k++ ) {
            j = Ngblist[k];
            //printf( "%i\n", P[j].Type );
            //
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
        Gprops[i].Head = Head[i];
        Gprops[i].Tail = Tail[i];
        Gprops[i].Len = Len[i];
    }

    myfree( Head );
    myfree( Tail );
    myfree( Len );
    myfree( Ngblist );

    writelog( "the maximum number of ngb: %i\n", ngbmax );
    mytimer_end();
    writelog( "FoF find groups ... done\n" );

}

void fof_compute_group_properties() {

    long p, i, k, p0;
    struct group_properties *g;
    double v2;

    writelog( "FoF compute groups properties ... \n" );
    mytimer_start();

    Ngroup = 0;
    for ( i=0; i<NumPart; i++ ) {
        if ( Gprops[i].Len >= All.FoFMinLen ) {
            Gprops[ Ngroup++ ] = Gprops[i];
        }
    }

    qsort( Gprops, Ngroup, sizeof( struct group_properties ), fof_compare_len );

    writelog( "The number of groups with at least %i particles: %i\n", All.FoFMinLen, Ngroup );
    writelog( "Largest group has %li particles\n", Gprops[0].Len );

    for ( i=0; i<Ngroup; i++ ) {

        g = &Gprops[i];

        p0 = p = g->Head;
        while( p>=0 ){

            for ( k=0, v2=0; k<3; k++ ){
                g->cm[k] += P[p].Mass *
                    ( PERIODIC_HALF( P[p].Pos[k]-P[p0].Pos[k] )
                      + P[p0].Pos[k] );
                g->vel[k] += P[p].Mass * P[p].Vel[k]; 
                v2 += SQR( P[p].Vel[k] );
            }

            g->v_mean += sqrt(v2);
            g->mass += P[p].Mass;
            g->npart[ P[p].Type ] ++;
            g->mass_table[ P[p].Type ] += P[p].Mass;

            p = FoFNext[p];

        }

        if ( g->mass == 0 ) {
            printf( "Gprops[%li] is zeros!!! Len: %li\n", i, g->Len );
            p = g->Head;
            while( p>=0 ){
                printf( "%g, %i\n", P[p].Mass, P[p].Type );
                p = FoFNext[p];
            }
            endrun(20181107);
        }

        for ( k=0; k<3; k++ ){
            g->cm[k] /= g->mass;
            g->vel[k] /= g->mass;
        }
        g->v_mean /= g->Len;

        p = g->Head;
        while( p>=0 ){
            for ( k=0, v2=0; k<3; k++ ) {
                vmax2( g->size[k], NGB_PERIODIC( P[p].Pos[k] - g->cm[k] ) );
                v2 += SQR( P[p].Vel[k] );
            }
            g->v_disp += SQR( sqrt(v2) - g->v_mean );

            for ( k=0, v2=0; k<3; k++ ) {
                v2 += SQR(  P[p].Vel[k] - g->vel[k] );
            }
            g->ek += 0.5 * P[p].Mass * v2 * Time3; 
            // to comoving energy,  v is comoving speed

            p = FoFNext[p];
        }
        g->vr200 = pow( g->mass /
            ( RhoCrit * 200 * 4.0 / 3.0 * PI ), 1.0/3.0 );

        g->v_disp /= g->Len;
        g->v_disp = sqrt( g->v_disp );

    }

    qsort( Gprops, Ngroup, sizeof( struct group_properties ), fof_compare_mass );

    mytimer_end();
    writelog( "FoF compute groups properties ... done\n" );

}

void check_fof( int gi, int flag ) {

    struct group_properties *g;
    long p;
    int i;
    for( i=0; i<100; i++)
        printf( "%li ", FoFNext[i] );
    printf( "\n" );
    g= &Gprops[gi];
    p = g->Head;
    while( p>=0 ) {
        for ( i=0; i<3; i++ )
            if ( PERIODIC_HALF( P[p].Pos[i] - g->cm[i] ) > g->size[i] ) {
                printf( "[flag: %i][%i]%g %g %g\n",
                        flag, i, P[p].Pos[i], g->cm[i], g->size[i]  );
                 endruns( "can't occur!"  );
            }
        p = FoFNext[p];
    }

}

#endif

void fof() {

#ifdef FOF
    long i, npart;
    int flag, num;
    double masstot, mass;
    char fn[ FILENAME_MAX ];


    if ( ThisTask_Local == 0 ) {
        writelog( "Start FoF ...\n" );
        sprintf( fn, "%s/fof_%03i.hdf5", All.FoFDir, SnapIndex );

        flag = 1;

        if ( access( fn, 0 ) != -1 ) {
            fof_read();
            flag = 0;
        }

        num = 0;

        MPI_Reduce( &flag, &num, 1, MPI_INT, MPI_SUM, 0, MpiComm_Master );
        MPI_Bcast( &num, 1, MPI_INT, 0, MpiComm_Master );

        create_dir( All.FoFDir );
        writelog( "%i snapshots Need to do FoF ...\n", num );
    }

    MPI_Bcast( &flag, 1, MPI_INT, 0, MpiComm_Local );

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


    if ( flag == 0 ) {

        if ( ThisTask_Local != 0 )
            mymalloc2( FoFNext, NumPart * sizeof( long ) );

        MPI_Bcast( &Ngroup, 1, MPI_INT, 0, MpiComm_Local );
        MPI_Bcast( FoFNext, NumPart*sizeof(long), MPI_BYTE, 0, MpiComm_Local );

        if ( ThisTask_Local != 0 )
            mymalloc2( Gprops, Ngroup * sizeof( struct group_properties ) );

        MPI_Bcast( Gprops, Ngroup*sizeof( struct group_properties ),
                MPI_BYTE, 0, MpiComm_Local );

        writelog( "%g\n", Gprops[0].mass);
        return;

    }

    mytimer_start();
    mymalloc2( Gprops, NumPart * sizeof( struct group_properties ) );
    mymalloc2( FoFNext, NumPart * sizeof( long ) );

    rhodm = (Omega0-All.OmegaBaryon) * 3 *  SQR(Hubble) / ( 8 * PI * G );
    mass = masstot / npart;
    LinkL = All.LinkLength * pow( mass / rhodm, 1.0/3 );
    //LinkL = 65.5483;
    /*
    printf( "%.10f %.10f %.10f %.10f %.20f %.10f %.10f %li \n",
         Omega0, All.OmegaBaryon,
         Hubble, G, rhodm, M_PI, masstot, npart );
         */

    writelog( "critical density of dark matter: %g\n"
        "comoving linking lenght %g\n", rhodm, LinkL );
    fof_find_groups();
    fof_compute_group_properties();

        writelog( "%g\n", Gprops[0].mass);
    if ( ThisTask_Local == 0 )
        fof_save();

    mytimer_end();
    writelog( "FoF ... done\n" );

#endif
}

void fof_free() {
#ifdef FOF
    writelog( "fof free\n" );
    myfree( Gprops );
    myfree( FoFNext );
#endif
}

void test_fof() {

#ifdef FOF_TEST
printf( "test_fof...\n" );
    tree_build();
    fof();

    check_fof( 0, 0 );

    endruns( "test_fof" );
#endif
}
