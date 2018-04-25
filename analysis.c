#include "allvars.h"

void init_analysis() {
    writelog( "initialize analysis...\n" );
    proj_k = All.ProjectDirection;
    proj_i = ( All.ProjectDirection + 1 ) % 3;
    proj_j = ( All.ProjectDirection + 2 ) % 3;
    slice();
    init_kernel_matrix();
    inte_ws = gsl_integration_workspace_alloc( GSL_INTE_WS_LEN );
    writelog( "initialize analysis... done.\n" );
    writelog( sep_str );
}

void free_analysis() {
    writelog( "free analysis ...\n" );
    gsl_integration_workspace_free( inte_ws );
    free_kernel_matrix();
    writelog( "free analysis ... done.\n" );
    writelog( sep_str );
}

void gas_density_slice() {
    int num, i;
    char buf[100];
    writelog( "gas density silce ...\n" );
    num = SliceEnd[0] - SliceStart[0];
    check_buffersize( (num + SQR(All.PicSize)) * sizeof( double ) );
    image.data = CommBuffer;
    image.img = (double *)CommBuffer + num;
    for ( i=SliceStart[0]; i<num; i++ ) {
        image.data[i] = SphP[i].Density;
    }

    create_dir( "./gas_density" );
    sprintf( buf, "./gas_density/%.2f.dat", All.RedShift );

    for ( i=0; i<SQR(All.PicSize); i++ ) {
        image.img[i] *= All.UnitMass_in_g / pow( All.UnitLength_in_cm, 2 );
    }
    image.DataMin *= All.UnitDensity_in_cgs;
    image.DataMax *= All.UnitDensity_in_cgs;
    image.GlobalDataMin *= All.UnitDensity_in_cgs;
    image.GlobalDataMax *= All.UnitDensity_in_cgs;

    image.ImgMin *= All.UnitMass_in_g / pow( All.UnitLength_in_cm, 2 );
    image.ImgMax *= All.UnitMass_in_g / pow( All.UnitLength_in_cm, 2 );
    image.GlobalImgMin *= All.UnitMass_in_g / pow( All.UnitLength_in_cm, 2 );
    image.GlobalImgMax *= All.UnitMass_in_g / pow( All.UnitLength_in_cm, 2 );

    make_slice_img( 0 );
    write_img( buf, 1 );

    writelog( "gas density silce ... done.\n" );
    writelog( sep_str );
}

int compare_gas_rho( const void *a, const void *b ){
    return ((((struct sph_particle_data* )a)->Density) < (((struct sph_particle_data*)b)->Density)) ? 1: -1;
}

void sort_gas_rho(){
    int i;
    qsort( (void*)SphP, N_Gas, sizeof( struct sph_particle_data ), &compare_gas_rho );
    for ( i=0; i<10; i++ ) {
        printf( "%g\n", SphP[i].Density * ( g2c.g / CUBE(g2c.cm) ) );
    }
    printf( "\n" );
    for ( i=0; i<10; i++ ) {
        printf( "%g\n", SphP[N_Gas-10+i].Density * ( g2c.g / CUBE(g2c.cm) ) );
    }
}

void vel_value() {
    writelog( "velocities value analysis ...\n" );
    FILE *fd;
    char buf[20];
    int i;
    double v;
    sprintf( buf, "vel_%i.txt", ThisTask );
    fd = fopen( buf, "w" );
    for ( i=0; i<N_Gas; i++ ) {
        v = sqrt( pow( P[i].Vel[0], 2 ) + pow( P[i].Vel[1], 2 ) + pow( P[i].Vel[2], 2 ) );
        if ( v > 1000 ) {
            fprintf( fd, "%i %g\n", P[i].ID, v );
        }
    }
    fclose( fd );
    writelog( "velocities value analysis ... done\n" );
    writelog( sep_str );
}

void test_id() {
    long id_max, id_min, i, num, offset;
    id_min = INT_MAX;
    id_max = INT_MIN;
    for ( i=0; i<N_Gas; i++ ) {
        id_max = ( P[i].ID > id_max ) ? P[i].ID : id_max;
        id_min = ( P[i].ID < id_min ) ? P[i].ID : id_min;
    }
    printf( "Gas MAX ID: %li, MIN ID: %li\n", id_max, id_min );
    offset = get_particle_offset( 4 );
    num = get_particle_offset( 4 );
    printf( "start offset: %li, num: %li\n", offset, num );
    if ( num == 0 ) return;
    id_min = INT_MAX;
    id_max = INT_MIN;
    for ( i=offset; i<offset+num; i++ ) {
        id_max = ( P[i].ID > id_max ) ? P[i].ID : id_max;
        id_min = ( P[i].ID < id_min ) ? P[i].ID : id_min;
    }
    printf( "Star MAX ID: %li, MIN ID: %li\n", id_max, id_min );
}

void calc_vl() {
    long i;
    double B;
    writelog( "compute Larmor frequency ...\n" );
    for ( i=0; i<N_Gas; i++ ) {
        B = sqrt( pow( SphP[i].B[0], 2 ) + pow( SphP[i].B[1], 2 ) + pow( SphP[i].B[2], 2 ) );
        SphP[i].vL = ELECTRON_CHARGE * B / ( 2 * PI * ELECTRON_MASS * LIGHT_SPEED );
    }
    writelog( "compute Larmor frequency ... done.\n" );
    writelog( sep_str );
}

void syn( double v ) {
    long i, xi, yi, PicSize;
    double B, Ub, vl, pmin, pmax, *img, dx, dy, d1, d2,
           com_dis, ang_dis, lum_dis, BoxSize, C, img_max, img_min, V,
           h, tmp, ang, beam;
    writelog( "compute synchrotron radiation ...\n" );
    pmin = DBL_MAX;
    pmax = DBL_MIN;
    PicSize = All.PicSize;
    BoxSize = All.BoxSize;
    for ( i=0; i<N_Gas; i++ ) {
        C = SphP[i].CRE_C0 * pow( SphP[i].Density, (All.Alpha-2)/3.0 );
        C = C * SphP[i].Density / ( ELECTRON_MASS /(g2c.g) );
        //printf( "%g\n", C );
        C /= CUBE( g2c.cm );
        //printf( "%g\n", C );
        B = sqrt( pow( SphP[i].B[0], 2 ) + pow( SphP[i].B[1], 2 ) + pow( SphP[i].B[2], 2 ) );
        Ub = B * B / ( 8 * PI );
        vl = SphP[i].vL;
        SphP[i].P = C * 2.0 / 3.0 * LIGHT_SPEED * Ub * THOMSON_CROSS_SECTION *
            pow( v / vl-1, (1-All.Alpha) / 2 ) / vl;
        pmin = ( SphP[i].P < pmin && SphP[i].P > 0 ) ? SphP[i].P : pmin;
        pmax = ( SphP[i].P > pmax ) ? SphP[i].P : pmax;
    }
    writelog( "pmin = %g, pmax = %g\n", pmin, pmax );

    check_buffersize( sizeof( double ) * SQR( PicSize ) );
    img = CommBuffer;
    memset( img, 0, sizeof( double ) * SQR( PicSize ) );

    com_dis = comoving_distance( header.time ) * ( g2c.cm );
    ang_dis = angular_distance( header.time ) * ( g2c.cm );
    lum_dis = luminosity_distance( header.time ) * ( g2c.cm );
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
        /*
        tmp = img[xi * PicSize + yi] += SphP[i].P * V / ( 4*PI*pow(lum_dis,2) ) / SQR(ang)
            * 1e25;
            */
        img_max = ( tmp > img_max ) ? tmp : img_max;
        img_min = ( tmp < img_min && tmp > 0 ) ? tmp : img_min;
    }
    writelog( "img_max = %g, img_min = %g\n", img_max, img_min );
}

void output_mag() {
    long i;
    FILE *fd;
    fd = fopen( "./B.txt", "w" );
    for ( i=0; i<N_Gas; i++ ) {
        fprintf( fd, "%g %g %g %g\n",
                P[i].Pos[0],
                P[i].Pos[1],
                P[i].Pos[2],
                log10(get_B( i ) ) );
    }
    fclose( fd );
}

void output_rho() {
    long i;
    FILE *fd;
    double rho_max, rho_min;
    fd = fopen( "./rho.txt", "w" );
    rho_max = -DBL_MAX;
    rho_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        fprintf( fd, "%g %g %g %g\n",
                P[i].Pos[0],
                P[i].Pos[1],
                P[i].Pos[2],
                SphP[i].Density );
        rho_max = ( SphP[i].Density > rho_max ) ? SphP[i].Density : rho_max;
        rho_min = ( SphP[i].Density < rho_min ) ? SphP[i].Density : rho_min;
    }
    printf( "rho_max = %g, rho_min = %g \n", rho_max, rho_min );
    fclose( fd );
}

void compute_temperature() {
    double yhelium, u, ne, mu, XH;
    int i;
    writelog( "compute gas temprature...\n" );
    XH = HYDROGEN_MASSFRAC;
    yhelium = ( 1 - XH ) / ( 4 * XH );
    for ( i=0; i<N_Gas; i++ ) {
        u = SphP[i].u * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
        ne = SphP[i].elec;
        mu = ( 1 + 4 * yhelium ) / ( 1 + yhelium + ne );
        SphP[i].Temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
    }
    writelog( "compute gas temprature... done.\n" );
    writelog( sep_str );
}

void gas_state() {
    double TempMin, TempMax, DensMin, DensMax,
           LogTempMin, LogTempMax, LogDensMin, LogDensMax,
           *img, DLogTemp, DLogDens, LogDens, LogTemp;
    int num, i, j, k, PicSize, N;
    char buf[100];
    TempMin = DensMin = DBL_MAX;
    TempMax = DensMax = DBL_MIN;
    writelog( "plot gas state...\n" );

    //PicSize_tmp = All.PicSize;
    PicSize = All.PicSize;

    check_buffersize( sizeof( double ) * N_Gas  + sizeof( double ) * SQR(PicSize) );

    img = CommBuffer;

    for ( i=0; i<N_Gas; i++ ) {
        TempMin = vmin( SphP[i].Temp, TempMin, 1 );
        TempMax = vmax( SphP[i].Temp, TempMax );
        DensMin = vmin( SphP[i].Density, DensMin, 1 );
        DensMax = vmax( SphP[i].Density, DensMax );
    }

    LogDensMin = log10( DensMin );
    LogDensMax = log10( DensMax );
    LogTempMin = log10( TempMin );
    LogTempMax = log10( TempMax );
    DLogTemp = ( LogTempMax - LogTempMin ) / PicSize;
    DLogDens = ( LogDensMax - LogDensMin ) / PicSize;
    writelog( "DensMin: %g, DensMax: %g, TempMin: %g, TempMax: %g\n"
            "LogDensMin: %g, LogDensMax: %g, LogTempMIn: %g, LogTempMax: %g\n",
            DensMin, DensMax, TempMin, TempMax,
            LogDensMin, LogDensMax, LogTempMin, LogTempMax );

    memset( img, 0, SQR( PicSize ) * sizeof( double ) );
    for ( k=0; k<N_Gas; k++ ) {
        LogTemp = SphP[k].Temp;
        //if ( LogTemp > 1e8 ) continue;
        LogTemp = ( LogDens == 0 ) ? LogTempMin : log10( LogTemp );
        LogDens = SphP[k].Density;
        LogDens = ( LogDens == 0 ) ? LogDensMin : log10( LogDens );

        i = ( LogTemp - LogTempMin ) / DLogTemp;
        i = check_picture_index( i );

        j = ( LogDens - LogDensMin ) / DLogDens;
        j = check_picture_index( j );

        img[ i*PicSize + j ]++;
 //       printf( "%g ", img[ i*PicSize +j ] ) ;
    }


    create_dir( "./gas_state" );
    sprintf( buf, "./gas_state/%.2f.dat", All.RedShift );

    memset( &image, 0, sizeof( struct image_struct ) );

    image.ImgMin = DBL_MAX;
    image.ImgMax = DBL_MIN;
    for ( i=0; i<SQR(PicSize); i++ ) {
        image.ImgMin = vmin( image.ImgMin, img[i], 1 );
        image.ImgMax = vmax( image.ImgMax, img[i] );
    }

    find_global_value( image.ImgMin, image.GlobalImgMin, MPI_DOUBLE, MPI_MIN );
    find_global_value( image.ImgMax, image.GlobalImgMax, MPI_DOUBLE, MPI_MAX );

    image.img = img;
    image.xmin = log10(DensMin * All.UnitDensity_in_cgs);
    image.xmax = log10(DensMax * All.UnitDensity_in_cgs);
    image.ymin = LogTempMin;
    image.ymax = LogTempMax;
    write_img( buf, 1 );

    //All.PicSize = PicSize_tmp;
    writelog( sep_str );
}

void gas_temperature_slice() {

    int num, i;
    char buf[100];
    writelog( "gas temperature silce ...\n" );
    num = SliceEnd[0] - SliceStart[0];
    check_buffersize( (num + SQR(All.PicSize)) * sizeof( double ) );
    image.data = CommBuffer;
    image.img = (double *)CommBuffer + num;
    for ( i=SliceStart[0]; i<num; i++ ) {
        image.data[i] = SphP[i].Temp;
    }

    create_dir( "./gas_temperature" );
    sprintf( buf, "./gas_temperature/%.2f.dat", All.RedShift );

    for ( i=0; i<SQR(All.PicSize); i++ ){
        image.img[i] *= All.UnitLength_in_cm;
    }

    image.ImgMin *=  All.UnitLength_in_cm;
    image.ImgMax *=  All.UnitLength_in_cm;
    image.GlobalImgMin *= All.UnitLength_in_cm;
    image.GlobalImgMax *= All.UnitLength_in_cm;

    make_slice_img( 0 );
    write_img( buf, 1 );

    writelog( "gas Temperature silce ... done.\n" );
    writelog( sep_str );

}

void analysis(){
    init_analysis();

    if ( All.GasTemperature ||
         All.GasState )
        compute_temperature();

    if ( All.GasState ) {
        gas_state();
    }

    if ( All.GasDensity )
        gas_density_slice();

    if ( All.GasTemperature )
        gas_temperature_slice();
    //tree_build( 1 );
    //tree_free();
    //fof( 1 );
    //fof_save_groups();
    //fof_free();
    //output_mag();
    //output_rho();
    //test_id();
    //vel_value();
    //sort_gas_rho();
    //calc_vl();
    //syn( 1.4e9 );
    //radio_halo();
    free_analysis();
}
