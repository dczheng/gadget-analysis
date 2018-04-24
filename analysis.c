#include "allvars.h"

FILE *fp_tmp;

void init_analysis() {
    print_log( "initialize analysis..." );
    slice();
    init_projection();
    init_kernel_matrix();
    init_plot();
    inte_ws = gsl_integration_workspace_alloc( GSL_INTE_WS_LEN );
    print_log( "initialize analysis... done." );
    print_log( sep_str );
}

void free_analysis() {
    print_log( "free analysis ..." );
    gsl_integration_workspace_free( inte_ws );
    free_kernel_matrix();
    free_plot();
    print_log( "free analysis ... done." );
    print_log( sep_str );
}

void magnetic_field_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];
    pli.h = All.SofteningTable[0];
    pli.log_flag = 1;
    sprintf( pli.data_name, "B" );
    sprintf( pli.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( pli.ylabel, "" );
    sprintf( pli.cb_label, "(10^x) G/Kpc^2");
    sprintf( pli.title, "" );
    pli.PicSize = All.PicSize;

    if ( BufferBytes < sizeof( double ) * num  + sizeof( double ) * SQR(pli.PicSize)) {
        printf( "BufferSize is too small !\n" );
        endrun( 20180424 );
    }

    pli.data = CommBuffer;
    pli.img = ( double* )CommBuffer + num;

    memset( pli.data, 0, sizeof(double) * num );
    pli.istart = SliceStart[0];
    pli.iend = SliceEnd[0];
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        pli.data[i-SliceStart[0]] = sqrt( pow( SphP[i].B[0], 2 ) +
                                                pow( SphP[i].B[1], 2 ) +
                                                pow( SphP[i].B[2], 2 ) );
    plot_slice();
}

void gas_rho_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    pli.h = All.SofteningTable[0];
    pli.log_flag = 1;
    sprintf( pli.data_name, "rho" );
    sprintf( pli.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( pli.ylabel, "" );
    sprintf( pli.title, "" );
    sprintf( pli.cb_label, "(10^x) gcm^{-3}Kpc^{-2}");
    pli.global_colorbar_flag = 1;
    pli.istart = SliceStart[0];
    pli.iend = SliceEnd[0];
    pli.PicSize = All.PicSize;

    if ( BufferBytes < sizeof( double ) * num  + sizeof( double ) * SQR(pli.PicSize)) {
        printf( "BufferSize is too small !\n" );
        endrun( 20180424 );
    }

    pli.data = CommBuffer;
    pli.img = ( double* )CommBuffer + num;

    memset( pli.data, 0, sizeof(double) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        pli.data[i-SliceStart[0]] = SphP[i].Density  *  (g2c.g / CUBE(g2c.cm) );
    plot_slice();
}

void mach_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    pli.h = All.SofteningTable[0];
    pli.log_flag = 1;
    sprintf( pli.data_name, "mn" );
    sprintf( pli.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( pli.ylabel, "" );
    sprintf( pli.title, "" );
    sprintf( pli.cb_label, "(10^x)");
    pli.istart = SliceStart[0];
    pli.iend = SliceEnd[0];
    pli.PicSize = All.PicSize;

    if ( BufferBytes < sizeof( double ) * num  + sizeof( double ) * SQR(pli.PicSize)) {
        printf( "BufferSize is too small !\n" );
        endrun( 20180424 );
    }

    pli.data = CommBuffer;
    pli.img = ( double* )CommBuffer + num;

    memset( pli.data, 0, sizeof(double) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        pli.data[i-SliceStart[0]] = SphP[i].MachNumber;
    plot_slice();
}

void divB_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    pli.h = All.SofteningTable[0];
    pli.log_flag = 1;
    sprintf( pli.data_name, "divB" );
    sprintf( pli.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( pli.ylabel, "" );
    sprintf( pli.title, "" );
    sprintf( pli.cb_label, "(10^x)");
    pli.istart = SliceStart[0];
    pli.iend = SliceEnd[0];
    pli.PicSize = All.PicSize;

    if ( BufferBytes < sizeof( double ) * num  + sizeof( double ) * SQR(pli.PicSize)) {
        printf( "BufferSize is too small !\n" );
        endrun( 20180424 );
    }

    pli.data = CommBuffer;
    pli.img = ( double* )CommBuffer + num;

    memset( pli.data, 0, sizeof(double) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        pli.data[i-SliceStart[0]] = SphP[i].divB;
    plot_slice();
}

void gas_vel_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    pli.h = All.SofteningTable[0];
    pli.log_flag = 1;
    sprintf( pli.data_name, "vel" );
    sprintf( pli.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( pli.ylabel, "" );
    sprintf( pli.title, "" );
    sprintf( pli.cb_label, "(10^x)");
    pli.istart = SliceStart[0];
    pli.iend = SliceEnd[0];
    pli.PicSize = All.PicSize;

    if ( BufferBytes < sizeof( double ) * num  + sizeof( double ) * SQR(pli.PicSize)) {
        printf( "BufferSize is too small !\n" );
        endrun( 20180424 );
    }

    pli.data = CommBuffer;
    pli.img = ( double* )CommBuffer + num;

    memset( pli.data, 0, sizeof(double) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        pli.data[i-SliceStart[0]] = sqrt( pow( P[i].Vel[0], 2 ) +
                                                pow( P[i].Vel[1], 2 ) +
                                                pow( P[i].Vel[2], 2 ) );
    plot_slice();
}

void hge_n_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    pli.h = All.SofteningTable[0];
    pli.log_flag = 1;
    pli.global_colorbar_flag = 1;
    sprintf( pli.data_name, "hge_n" );
    sprintf( pli.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( pli.ylabel, "" );
    sprintf( pli.title, "" );
    sprintf( pli.cb_label, "(10^x)cm^{-3}Kpc^{-2}");
    pli.istart = SliceStart[0];
    pli.iend = SliceEnd[0];
    pli.PicSize = All.PicSize;

    if ( BufferBytes < sizeof( double ) * num  + sizeof( double ) * SQR(pli.PicSize)) {
        printf( "BufferSize is too small !\n" );
        endrun( 20180424 );
    }

    pli.data = CommBuffer;
    pli.img = ( double* )CommBuffer + num;

    memset( pli.data, 0, sizeof(double) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        pli.data[i-SliceStart[0]] = SphP[i].CRE_n0 * SphP[i].Density * ( g2c.g / CUBE(g2c.cm) ) / ELECTRON_MASS;
    plot_slice();
}

void cr_n_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    pli.h = All.SofteningTable[0];
    pli.log_flag = 1;
    pli.global_colorbar_flag = 1;
    sprintf( pli.data_name, "cr_n" );
    sprintf( pli.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( pli.ylabel, "" );
    sprintf( pli.title, "" );
    sprintf( pli.cb_label, "(10^x)cm{-3}Kpc^{-2}");
    pli.istart = SliceStart[0];
    pli.iend = SliceEnd[0];
    pli.PicSize = All.PicSize;

    if ( BufferBytes < sizeof( double ) * num  + sizeof( double ) * SQR(pli.PicSize)) {
        printf( "BufferSize is too small !\n" );
        endrun( 20180424 );
    }

    pli.data = CommBuffer;
    pli.img = ( double* )CommBuffer + num;

    memset( pli.data, 0, sizeof(double) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        pli.data[i-SliceStart[0]] = SphP[i].CR_n0 * SphP[i].Density * ( g2c.g / CUBE(g2c.cm) ) / PROTONMASS;
    plot_slice();
}

void cr_e_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    pli.h = All.SofteningTable[0];
    pli.log_flag = 1;
    pli.global_colorbar_flag = 0;
    sprintf( pli.data_name, "cr_e" );
    sprintf( pli.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( pli.ylabel, "" );
    sprintf( pli.title, "" );
    sprintf( pli.cb_label, "(10^x)erg cm^{-3} Kpc^{-2}");
    pli.istart = SliceStart[0];
    pli.iend = SliceEnd[0];
    pli.PicSize = All.PicSize;

    if ( BufferBytes < sizeof( double ) * num  + sizeof( double ) * SQR(pli.PicSize)) {
        printf( "BufferSize is too small !\n" );
        endrun( 20180424 );
    }

    pli.data = CommBuffer;
    pli.img = ( double* )CommBuffer + num;

    memset( pli.data, 0, sizeof(double) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        pli.data[i-SliceStart[0]] = SphP[i].CR_E0 * SphP[i].Density  * (g2c.erg / CUBE(g2c.cm));
    plot_slice();
}

void hge_e_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    pli.h = All.SofteningTable[0];
    pli.log_flag = 1;
    pli.global_colorbar_flag = 0;
    sprintf( pli.data_name, "hge_e" );
    sprintf( pli.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( pli.ylabel, "" );
    sprintf( pli.title, "" );
    sprintf( pli.cb_label, "(10^x)erg cm^{-3} Kpc^{-2}");
    pli.istart = SliceStart[0];
    pli.iend = SliceEnd[0];
    pli.PicSize = All.PicSize;

    if ( BufferBytes < sizeof( double ) * num  + sizeof( double ) * SQR(pli.PicSize)) {
        printf( "BufferSize is too small !\n" );
        endrun( 20180424 );
    }

    pli.data = CommBuffer;
    pli.img = ( double* )CommBuffer + num;

    memset( pli.data, 0, sizeof(double) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        pli.data[i-SliceStart[0]] = SphP[i].CRE_E0 * SphP[i].Density * (g2c.erg / CUBE(g2c.cm) );
    plot_slice();
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
    print_log( "velocities value analysis ..." );
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
    print_log( "velocities value analysis ... done" );
    print_log( sep_str );
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
    print_log( "compute Larmor frequency ..." );
    for ( i=0; i<N_Gas; i++ ) {
        B = sqrt( pow( SphP[i].B[0], 2 ) + pow( SphP[i].B[1], 2 ) + pow( SphP[i].B[2], 2 ) );
        SphP[i].vL = ELECTRON_CHARGE * B / ( 2 * PI * ELECTRON_MASS * LIGHT_SPEED );
    }
    print_log( "compute Larmor frequency ... done." );
    print_log( sep_str );
}

void vl_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    pli.h = All.SofteningTable[0];
    pli.log_flag = 1;
    pli.global_colorbar_flag = 0;
    sprintf( pli.data_name, "vl" );
    sprintf( pli.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( pli.ylabel, "" );
    sprintf( pli.title, "" );
    sprintf( pli.cb_label, "(10^x)(MHz)");
    pli.istart = SliceStart[0];
    pli.iend = SliceEnd[0];
    pli.PicSize = All.PicSize;

    if ( BufferBytes < sizeof( double ) * num  + sizeof( double ) * SQR(pli.PicSize)) {
        printf( "BufferSize is too small !\n" );
        endrun( 20180424 );
    }

    pli.data = CommBuffer;
    pli.img = ( double* )CommBuffer + num;

    memset( pli.data, 0, sizeof(double) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        pli.data[i-SliceStart[0]] = SphP[i].vL / 1e6;
    plot_slice();
}

void syn( double v ) {
    long i, xi, yi, PicSize;
    double B, Ub, vl, pmin, pmax, *img, dx, dy, d1, d2,
           com_dis, ang_dis, lum_dis, BoxSize, C, img_max, img_min, V,
           h, tmp, ang, beam;
    print_log( "compute synchrotron radiation ..." );
    pmin = DBL_MAX;
    pmax = DBL_MIN;
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
    sprintf( LogBuf, "pmin = %g, pmax = %g", pmin, pmax );
    print_log( LogBuf );

    if ( BufferBytes < sizeof( double ) * SQR( All.PicSize ) ) {
        printf( "BufferSize is too small !\n" );
        endrun( 20180424 );
    }
    img = CommBuffer;
    memset( img, 0, sizeof( double ) * SQR( All.PicSize ) );

    com_dis = comoving_distance( header.time ) * ( g2c.cm );
    ang_dis = angular_distance( header.time ) * ( g2c.cm );
    lum_dis = luminosity_distance( header.time ) * ( g2c.cm );
    PicSize = All.PicSize;
    BoxSize = All.BoxSize;
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
    sprintf( LogBuf, "img_max = %g, img_min = %g\n", img_max, img_min );
    print_log( LogBuf );
    pli.img = img;
    sprintf( pli.data_name, "syn" );
    pli.log_flag = 1;
    sprintf( pli.xlabel, "%.2g deg", BoxSize * (g2c.cm) / ang_dis / PI * 180);
    sprintf( pli.ylabel, "%.2g deg", BoxSize * (g2c.cm) / ang_dis / PI * 180);
    sprintf( pli.title, "1.4GHz(z=%.2f)", All.RedShift );
    //sprintf( pli.cb_label, "10^x(mJy arcmin^{-2})" );
    sprintf( pli.cb_label, "10^x(mJy beam^{-1})" );
    plot_imshow();
    print_log( "compute synchrotron radiation ... done." );
    print_log( sep_str );
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

void gas_state() {
    double yhelium, mu, XH=HYDROGEN_MASSFRAC, T,
           Tmin, Tmax, rhomin, rhomax, *img;
    int num, i, k, PicSize;
    Tmin = rhomin = DBL_MAX;
    Tmax = rhomax = DBL_MIN;

    yhelium = ( 1 - XH ) / ( 4 * XH );

    img = malloc( SQR(PicSize) * sizeof( double ) );
    for ( i=0; i<N_Gas; i++ ) {
        mu = ( 1 + 4 * yhelium ) / ( 1 + yhelium + SphP[i].elec );
        T = GAMMA_MINUS1 / BOLTZMANN * SphP[i].u * PROTONMASS * mu;
        Tmin = ( T < Tmin ) ? T : Tmin;
        Tmax = ( T > Tmax ) ? T : Tmax;
        rhomin = ( SphP[i].Density < rhomin ) ? SphP[i].Density : rhomin;
        rhomax = ( SphP[i].Density > rhomax ) ? SphP[i].Density : rhomax;
    }

    printf( "rhomin: %g, rhomax: %g, Tmin: %g, Tmax: %g\n",
            rhomin, rhomax, Tmin, Tmax );

    rhomin = log10( rhomin );
    rhomax = log10( rhomax );
    Tmin = log10( Tmin );
    Tmax = log10( Tmax );
}


void analysis(){
    init_analysis();
    //gas_state();
    //tree_build( 1 );
    //tree_free();
    //fof( 1 );
    //fof_save_groups();
    //fof_free();
    gas_rho_slice();
    //magnetic_field_slice();
    //output_mag();
    //output_rho();
    //test_id();
    //divB_slice();
    //gas_vel_slice();
    //mach_slice();
    //hge_n_slice();
    //cr_n_slice();
    //hge_e_slice();
    //cr_e_slice();
    //vel_value();
    //sort_gas_rho();
    //calc_vl();
    //syn( 1.4e9 );
    //radio_halo();
    free_analysis();
}
