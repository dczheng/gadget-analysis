#include "allvars.h"

#ifndef ALTRAD
#ifdef RADSPEC

void d2PdVdv_qmax() {

    SphParticleData part;
    double nu_min, nu_max, nu, dlognu, *J, *J0, logqmax_min, logqmax_max, dlogqmax;
    int i, N, qmaxn, k;
    FILE *fd;

    part.CRE_C = CUBE(g2c.cm) * ( guc.m_e ) / RhoBaryon;
    part.CRE_Alpha = 2.01;
    part.CRE_qmin = 1;
    part.CRE_qmax = 1e8;
    part.Density = RhoBaryon;

    part.B[0] = 1e-6;
    //part.B[0] = sqrt( SQR(B) - SQR(BCMB0) * pow( Time, -2 ) );
    part.B[1] = part.B[2] = 0;

    N = 100;
    nu_min = 1e6;
    nu_max = 1e9;
    dlognu = log10( nu_max/nu_min ) / ( N-1 );

    qmaxn = 10;
    logqmax_min = 4;
    logqmax_max = 5;
    dlogqmax = (logqmax_max - logqmax_min) / ( qmaxn-1 );

    mymalloc2( J, N * sizeof(double) );

    if ( ThisTask == 0 ) {
        mymalloc2( J0, N * sizeof(double) );
        fd = myfopen( "w", "./d2PdVdv_qmax.dat");
        fprintf( fd, "0 " );

        for( i=0; i<N; i++ ) {
            nu = log10(nu_min) + i * dlognu;
            nu = pow( 10, nu );
            fprintf( fd, "%g ", nu );
        }

        fprintf( fd, "\n" );
    }

    for( k=0; k<qmaxn; k++ ) {
        part.CRE_qmax = pow( 10, logqmax_min + k*dlogqmax );

        //part.CRE_qmax = qmax_max;

        if ( ThisTask == 0 ) {
            fprintf( fd, "%g ",
                part.CRE_qmax );
            printf( "qmax: %g \n",
                part.CRE_qmax );
        }

        /*
        printf(  "%g %g %g %g %g\n",
                part.CRE_C /( CUBE(g2c.cm) * ( cuc.m_e/(g2c.g) )),
                part.CRE_Alpha,
                part.CRE_qmin,
                part.CRE_qmax,
                part.B[0] );
                */

        for( i=0; i<N; i++ ) {
            if ( i % NTask != ThisTask )
                continue;

           // printf( "Task %i, i: %i\n", ThisTask_Local, i );

            nu = log10(nu_min) + i * dlognu;
            nu = pow( 10, nu );
            J[i] = particle_radio2( nu, &part ) / ( 4.0/3.0 * PI * CUBE( SofteningTable[0] * g2c.cm ) );
            //printf( "nu: %g, P: %g\n", nu, P );
            //break;
        }

        MPI_Reduce( J, J0, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

        if ( ThisTask == 0 ) {
            for( i=0; i<N; i++ )
                fprintf( fd, "%g ", J0[i] );
            fprintf( fd, "\n" );
            fflush( fd );
        }

       // break;

    }

    if ( ThisTask == 0 ) {
        myfree( J0 );
        fclose( fd );
    }
    myfree( J );
    do_sync( "" );

    writelog( "d2PdVdv_qmax done.\n" );

}

double particle_df( double p, void *params ) {

    double *pa;
    pa = params;

    if ( p < pa[2] || p > pa[3] )
        return 0;

    return pa[0] * pow( p, -pa[1] );

}

double particle_radio2( double nu,  SphParticleData *part ) {

    double r, params[4], B;

    params[0] = part->CRE_C;
    if ( part->CRE_qmin <= 0 )
        return 0;
    /*
    params[0] = part->CRE_n * ( part->CRE_Alpha - 1 ) /
        pow( part->CRE_qmin, 1-part->CRE_Alpha );
        */
    params[0] = params[0] * part->Density / guc.m_e;
    params[0] /= CUBE( g2c.cm );
    params[0] /= Time3;
    /*
    printf( "n: %g, a: %g, qmin: %g, %g\n",
            part->CRE_n, part->CRE_Alpha, part->CRE_qmin,
            params[0] );
            */

    params[1] = part->CRE_Alpha;
    params[2] = part->CRE_qmin;
    params[3] = part->CRE_qmax;

    /*
    printf( "c: %g, a: %g, qmin: %g, qmax: %g\n", params[0],
            params[1], params[2], params[3] );
            */

    if (  params[ 0] *
            params[1] *
            params[2] *
            params[3] == 0 )
        return 0;

    B = SQR(part->B[0]) + SQR(part->B[1]) + SQR(part->B[2]);
    B = sqrt( B );
    //if ( B > 1e-4 )
    //    B = 1e-4;
    r = radio( &particle_df, params, B, nu, params[2], params[3], 1e-2 );
    return r;

}

double particle_radio( double nu, long i ) {

    return particle_radio2( nu, &SphP[i] );

}

void save_particle_radio() {

    char fn[ FILENAME_MAX ];
    int ndims;
    hsize_t dims[2];

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

void compute_particle_radio() {

    double numin, numax, dlognu, nu, tt, V;
    int nuN, j, num, flag;
    long i, t;
    char fn[ FILENAME_MAX ];

    writelog( "Start compute particle radio ... \n" )

    init_compute_F();

    mymalloc1_shared( PartRad, sizeof(double)*N_Gas*All.NuNum, sizeof(double), MpiWin_PartRad );

    if ( ThisTask_Local == 0 ) {
        //printf( "Task: %i enter %s\n", ThisTask, __FUNCTION__ );

        sprintf( fn, "%s/rad_%03i.hdf5", All.RadDir, SnapIndex );
        flag = 1;

        if ( access( fn, 0 ) != -1 ) {
            if ( read_particle_radio() )
                flag = 0;
        }

        MPI_Reduce( &flag, &num, 1, MPI_INT, MPI_SUM, 0, MpiComm_Master );
        MPI_Bcast(  &num, 1, MPI_INT, 0, MpiComm_Master );

        if ( num == NTask_Master )
            create_dir( All.RadDir );

        writelog( "%i Task Need to compute radio.\n", num );

    }

    MPI_Bcast(  &num, 1, MPI_INT, 0, MpiComm_Local );
    MPI_Bcast(  &flag, 1, MPI_INT, 0, MpiComm_Local );

    if ( num == 0 )
        return;

#ifdef TABF
        init_tab_F();
#endif

    if ( flag == 0 ) {
#ifdef TABF
        free_tab_F();
#endif
        return;
    }


    nuN = All.NuNum;
    numin = All.NuMin;
    numax = All.NuMax;

    dlognu = log( numax/numin ) / ( nuN - 1 );

    t = N_Gas / 10;
    for( i=0; i<N_Gas; i++ ) {

        if ( i % t == 0 )
            writelog( "[%10li] [%10li] [%5.1f%%]\n", i, N_Gas, ((double)(i)) / N_Gas * 100 );

        if ( i % NTask_Local != ThisTask_Local )
            continue;

        flag = 0;
        V = P[i].Mass / SphP[i].Density * CUBE( g2c.cm * Time );
        //V = SphP[i].Hsml * CUBE( g2c.cm * Time );
        for( j=0; j<nuN; j++ ) {

            nu = exp( log(numin) + j * dlognu ) * 1e6;

            //printf( "nu: %g\n", nu );

            tt = PartRad[ i*All.NuNum + j ] = particle_radio( nu/Time, i ) * V;

            //tt  = tt * ( 4.0/3.0 * PI * CUBE(SofteningTable[0]*g2c.cm*Time) );
            /*
            tt = tt / ( 4.0 * PI * SQR( LumDis*g2c.cm ) ) * 1e26; // to mJy
            
            if ( tt < 1e-10 ) {
                flag = 1;    // ignore particle with very weak radio emission.
                break;
            }
            */

        }

        if ( flag )
            for( j=0; j<nuN; j++ )
            PartRad[ i*All.NuNum + j ] = 0;

    }


#ifdef TABF
        free_tab_F();
#endif

    if ( ThisTask_Local == 0 )
        save_particle_radio();

}

void free_particle_radio() {

    myfree_shared( MpiWin_PartRad );

}

void output_radio_inte() {

    int i, N;
    double logqmax, logqmin, dlogq, *q, *r_local, *r, p[4];
    struct radio_inte_struct ri;
    FILE *fd;

    logqmin = -1;
    logqmax = 10;
    N = 100;

    dlogq = ( logqmax - logqmin ) / ( N-1 );

    mymalloc2( q, N * sizeof( double ) );
    mymalloc2( r_local, N * sizeof( double ) );

    if ( ThisTask_Local == 0 ) {
        mymalloc2( r, N * sizeof( double ) );
    }

    for( i=0; i<N; i++ )
        q[i] = pow(10, logqmin + i * dlogq);

    p[0] = 1 * CUBE(g2c.cm) * ( cuc.m_e/(g2c.g) );
    p[1] = 2.01;
    p[2] = 1;
    p[3] = 1e4;

    ri.f = particle_df;
    ri.params = &p;
    ri.B = 1e-6;
    ri.nu = 1e4;

    for( i=0; i<N; i++ ) {
        if ( i % NTask_Local != ThisTask_Local )
            continue;
        r_local[i] = radio_inte( q[i], &ri );
    }

    MPI_Reduce( r_local, r, N, MPI_DOUBLE, MPI_SUM, 0, MpiComm_Local );

    if ( ThisTask_Local == 0 ) {

        fd = myfopen( "w", "radio_inte.dat");

        for( i=0; i<N; i++ )
            fprintf( fd, "%g %g\n", q[i], r[i] );

        fclose( fd );

        myfree( r );

    }

    myfree( r_local );
    myfree( q );
}
#endif
#endif

void test_part_radio() {


//#define TEST_PART_RAIO
#ifdef TEST_PART_RAIO
    HubbleParam = 0.68;
    Redshift = 0.01;
    Time = (1/(Redshift+1));
    Omega0 = 0.302;
    OmegaLambda = 0.698;
    OmegaBaryon = 0.0471;
    set_units();
    compute_cosmo_quantities();

#ifndef ALTRAD

    SphParticleData part;
    double numin, numax, nu, p, dlognu, dlogp, dfp, dfnu;
    int N, i;
    FILE *fd;

    part.Density = 1;
    //part.CRE_C = 7.316051e-13;
    part.CRE_C = 1;
    //part.CRE_qmin = 2.483418e+02;
    part.CRE_qmin = 0.1;
    //part.CRE_qmax = 1.520993e+06;
    part.CRE_qmax = 1e+06;
    part.CRE_Alpha = 2.1;
    part.B[0] = part.B[1] = 0;
    //part.B[2] = 5.504009e-07;
    part.B[2] = 1e-7;
    N = 100;
    numin = 1e7;
    numax = 1.5e9;

    init_compute_F();

#ifdef TABF
    init_tab_F();
#endif

    if ( ThisTask == 0 ) {
        dlogp = log( part.CRE_qmax * 10 / (part.CRE_qmin / 10) ) / N;
        dlognu = log( numax / numin ) / N;
    
        fd = fopen( "./test_part_radio.dat", "w" );
        for ( i=0; i<N; i++ ) {
            p = part.CRE_qmin / 10  * exp( i * dlogp );
            nu = numin * exp( i * dlognu );
            dfp = CRE_F( part.CRE_C, part.CRE_Alpha, part.CRE_qmin,
                 part.CRE_qmax, p );
            dfnu = particle_radio2( nu,  &part );
            fprintf( fd, "%e %e %e %e\n", p, dfp, nu, dfnu );
        }

        fclose( fd );
    }

#ifdef TABF
free_tab_F();
#endif
#else

    double vL, B;
    B = 1e-6;
    vL = B / ( 2 * PI ) * cuc.e_mec;

if ( ThisTask == 0 )
    printf( "vL: %g\n", vL );

#endif
    
    do_sync( "" );
    endruns("test_part_radio");
#else
    return;
#endif

}

double get_particle_radio( long p, int i ) {

#ifndef ALTRAD
    return PartRad[ p*All.NuNum + i ];
#else
    double vL, B, a, pmin, pmax, C, t, r, V, nu;
    B = get_B( p );
    vL = B / ( 2*PI ) * cuc.e_mec;
    C = SphP[p].CRE_C * SphP[p].Density / (guc.m_e * CUBE(g2c.cm) * Time3);
    a = (SphP[p].CRE_Alpha-1)/2.0;
    pmin = SphP[p].CRE_qmin;
    pmax = SphP[p].CRE_qmax;
    V = P[p].Mass / SphP[p].Density * CUBE( g2c.cm * Time );
    nu = All.NuMin * exp( i * ( log(All.NuMax/All.NuMin) / All.NuNum ) ) * 1e6;

    t = sqrt(nu/vL);
    if ( t<pmin || t>pmax )
        return 0;

    r = 2.0/3.0 * cuc.sigma_t * C * cuc.c * SQR(B) / ( 8*PI * vL ) * pow( nu/vL, -a ); 
    //printf( "%g %g %g %g %g %g\n", nu, vL, t, pmin, pmax, r  );
    r *= V;
    return r;

#endif
}
