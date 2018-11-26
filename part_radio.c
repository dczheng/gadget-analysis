#include "allvars.h"

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
    params[0] = params[0] * part->Density / ( cuc.m_e /(g2c.g) );
    params[0] /= CUBE( g2c.cm );

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

    B = pow( part->B[0], 2 ) + pow( part->B[1], 2 ) + pow( part->B[2], 2 );
    B += SQR( BCMB0 ) * pow( All.Time, -2 );
    B = sqrt( B );

    //printf( "B: %g\n", B );

    r = radio( &particle_df, params, B, nu, params[2], params[3] );

    r = r * ( 4.0/3.0 * PI * CUBE( All.SofteningTable[0] * g2c.cm ) );

    return r;

}

double particle_radio( double nu, long i ) {

    return particle_radio2( nu, &SphP[i] );

}

void save_particle_radio() {

    char fn[ FILENAME_MAX ];
    double *buf;
    int ndims, j;
    long i;
    hsize_t dims[2];

    hid_t hdf5_file, hdf5_dataset, hdf5_dataspace, hdf5_attribute, hdf5_type;

    sprintf( fn, "%s/rad_%.2f.hdf5", All.RadDir, All.RedShift );

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

    mymalloc1( buf, sizeof(double) * N_Gas * All.NuNum );

    for ( i=0; i<N_Gas; i++ )
        for ( j=0; j<All.NuNum; j++ )
            buf[ i * All.NuNum + j ] = SphP[i].Rad[j];

    hdf5_dataspace = H5Screate_simple( ndims, dims, NULL );
    hdf5_type = H5Tcopy( H5T_NATIVE_DOUBLE );
    hdf5_dataset = H5Dcreate( hdf5_file, "Radio", hdf5_type, hdf5_dataspace, H5P_DEFAULT );
    H5Dwrite( hdf5_dataset, hdf5_type, hdf5_dataspace, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    H5Tclose( hdf5_type );
    H5Sclose( hdf5_dataspace );

    myfree( buf );

    H5Fclose( hdf5_file );


}

int read_particle_radio() {

    double *buf, NuMin, NuMax;
    int nuN, j;
    long i;
    char fn[ FILENAME_MAX ];
    hid_t hdf5_file, hdf5_dataset, hdf5_attribute, hdf5_type;

    writelog( "read radio ...\n" );

    sprintf( fn, "%s/rad_%.2f.hdf5", All.RadDir, All.RedShift );

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

    mymalloc1( buf, sizeof(double) * N_Gas * All.NuNum );

    hdf5_type = H5Tcopy( H5T_NATIVE_DOUBLE );
    hdf5_dataset = H5Dopen( hdf5_file, "Radio" );
    H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );
    H5Dclose( hdf5_dataset );
    H5Tclose( hdf5_type );

    for( i=0; i<N_Gas; i++ ) {

        SphP[i].Rad = malloc( sizeof( double ) * All.NuNum );
        for( j=0; j<All.NuNum; j++ )
            SphP[i].Rad[j] = buf[ i * All.NuNum +j ];

    }

    myfree( buf );

    H5Fclose( hdf5_file );

    return 1;

}

void compute_particle_radio() {

    double numin, numax, dlognu, nu;
    int nuN, signal, j, num, flag;
    long i;
    char fn[ FILENAME_MAX ];

    writelog( "Start compute particle radio ... \n" )

    sprintf( fn, "%s/rad_%.2f.hdf5", All.RadDir, All.RedShift );

    flag = 1;

    if ( access( fn, 0 ) != -1 ) {

        if ( read_particle_radio() )
            flag = 0;

    }

    MPI_Reduce( &flag, &num, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Bcast(  &num, 1, MPI_INT, 0, MPI_COMM_WORLD );


    if ( num == NTask )
        create_dir( All.RadDir );

    writelog( "%i Task Need to compute radio.\n", num );

    if ( flag == 0 )
        return;

#ifdef RADIO_F_INTERP
        init_tab_F();
#endif

    nuN = All.NuNum;
    numin = All.NuMin;
    numax = All.NuMax;

    dlognu = log( numax/numin ) / ( nuN - 1 );

    signal = N_Gas / 10;

    for( i=0; i<N_Gas; i++ ) {


        if ( i % signal == 0 )
            writelog( "[%10li] [%10li] [%5.1f%%]\n", i, N_Gas, ((double)(i)) / N_Gas * 100 );

        SphP[i].Rad = malloc( sizeof( double ) * nuN );

        for( j=0; j<nuN; j++ ) {

            nu = exp( log(numin) + j * dlognu ) * 1e6;

            //printf( "nu: %g\n", nu );

            SphP[i].Rad[j] = particle_radio( nu, i );

        }

    }


#ifdef RADIO_F_INTERP
        free_tab_F();
#endif

    save_particle_radio();

    //endrun( 20181027 );
}

void free_particle_radio() {

    long i;
    for( i=0; i<N_Gas; i++ )
        free( SphP[i].Rad );

}

void test_qmax() {

    SphParticleData part;
    double nu_min, nu_max, nu, dlognu, P, qmax_min, qmax_max, dlogqmax, B;
    int i, N, qmaxn, k;
    FILE *fd;

    part.CRE_C = 1 * CUBE(g2c.cm) * ( cuc.m_e/(g2c.g) );
    part.CRE_Alpha = 2.01;
    part.CRE_qmin = 1;
    part.CRE_qmax = 1e5;
    part.Density = 1;
    B = BCMB0;
    part.B[0] = sqrt( SQR(B) - SQR(BCMB0) * pow( All.Time, -2 ) );
    part.B[1] = part.B[2] = 0;

    N = 100;
    nu_min = 1e6;
    nu_max = 1e9;
    dlognu = log10( nu_max/nu_min ) / ( N-1 );

    qmaxn = 5;
    qmax_min = 5e3;
    qmax_max = 1e5;
    dlogqmax = log10( qmax_max/qmax_min ) / ( qmaxn-1 );

    fd = fopen( "./qmax_test.dat", "w" );

    fprintf( fd, "0 " );

    for( i=0; i<N; i++ ) {
        nu = log10(nu_min) + i * dlognu;
        nu = pow( 10, nu );
        fprintf( fd, "%g ", nu );
    }

    fprintf( fd, "\n" );

    for( k=0; k<qmaxn; k++ ) {
        part.CRE_qmax = log10( qmax_min ) + k*dlogqmax;
        part.CRE_qmax = pow( 10, part.CRE_qmax );

        //part.CRE_qmax = qmax_max;

        fprintf( fd, "%g ",
                part.CRE_qmax );

        printf(  "%g %g %g %g %g\n",
                part.CRE_C /( CUBE(g2c.cm) * ( cuc.m_e/(g2c.g) )),
                part.CRE_Alpha,
                part.CRE_qmin,
                part.CRE_qmax,
                part.B[0] );

        for( i=0; i<N; i++ ) {
            nu = log10(nu_min) + i * dlognu;
            nu = pow( 10, nu );
            P = particle_radio2( nu, &part );
            //printf( "nu: %g, P: %g\n", nu, P );
            //break;
            fprintf( fd, "%g ", P );
        }

        fprintf( fd, "\n" );

        //break;

    }

    fclose( fd );

}

void test_radio() {


    All.Time = 1;
    All.HubbleParam = 0.7;
    set_units();
    put_block_line;

#ifdef RADIO_F_INTERP
    init_tab_F();
#endif



    if ( ThisTask == 0 ) {
        test_F();
        //test_qmax();
    }

    MPI_Barrier( MPI_COMM_WORLD );
#ifdef RADIO_F_INTERP
    free_tab_F();
#endif
    endrun(20181004);

}

/*
void compute_radio( double v ) {

    mymalloc( img, sizeof( double ) * SQR( PicSize ) );
    memset( img, 0, sizeof( double ) * SQR( PicSize ) );

    dy = dx = BoxSize / PicSize;
    h = All.SofteningTable[0] * ( g2c.cm );
    V = 4.0 / 3.0 * PI * pow( h, 3 );
    ang = h / ang_dis / PI * 180.0 * 60;
    beam = pow( 10.0 / 60, 2.0 ) / ( 4.0 * log(2) );
    img_max = DBL_MIN;
    img_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        xi = (int)( P[i].Pos[0] / dx );
        yi = (int)( P[i].Pos[1] / dy );
        tmp = img[xi * PicSize + yi] += SphP[i].P * V / ( 4*PI*pow(lum_dis,2) ) / SQR(ang)
            * beam * 1e25;
        //tmp = img[xi * PicSize + yi] += SphP[i].P * V / ( 4*PI*pow(lum_dis,2) ) / SQR(ang)
        //    * 1e25;
        img_max = ( tmp > img_max ) ? tmp : img_max;
        img_min = ( tmp < img_min && tmp > 0 ) ? tmp : img_min;
    }
}
*/

/*
double particle_radio( double nu, long i ) {

    C = SphP[i].CRE_C;
    //printf( "%g\n", C );
    C = C * SphP[i].Density / ( cuc.m_e /(g2c.g) );
    C /= CUBE( g2c.cm );

    B = sqrt( pow( SphP[i].B[0], 2 ) + pow( SphP[i].B[1], 2 ) + pow( SphP[i].B[2], 2 ) );
    Ub = B * B / ( 8 * PI );
    Ub += SQR( BCMB0 ) * pow( All.Time, -4 ) / ( 8*PI );

    nuL = cuc.e * B / ( 2 * PI * cuc.m_e * cuc.c );

    if ( sqrt( nu/nuL-1 ) < SphP[i].CRE_qmin ||
         sqrt( nu/nuL-1 ) > SphP[i].CRE_qmax )
        return 0;

    P = C * 2.0 / 3.0 * cuc.c * Ub * THOMSON_CROSS_SECTION *
        pow( nu / nuL-1, (1-SphP[i].CRE_Alpha) / 2 ) / nuL;

    //printf( "%g\n", P );
    //
    P = P * ( 4.0/3.0 * PI * CUBE( All.SofteningTable[0] * g2c.cm ) );

    // Unit: erg / Hz / s
    //
    //printf( "Ub: %g, nuL: %g, P: %g\n",  Ub, nuL, P );

    printf( "P: %g, P2: %g\n", P, P2 );

    endrun( 20181005 );

    return P;

}
*/




