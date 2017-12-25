#include "allvars.h"

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
    plot_info.h = All.SofteningTable[0];
    plot_info.log_flag = 1;
    sprintf( plot_info.data_name, "B" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( plot_info.ylabel, "" );
    sprintf( plot_info.cb_label, "(10^x) G/Kpc^2");
    sprintf( plot_info.title, "" );
    plot_info.data = malloc( sizeof(double) * num );
    memset( plot_info.data, 0, sizeof(double) * num );
    plot_info.istart = SliceStart[0];
    plot_info.iend = SliceEnd[0];
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        plot_info.data[i-SliceStart[0]] = sqrt( pow( SphP[i].B[0], 2 ) +
                                                pow( SphP[i].B[1], 2 ) +
                                                pow( SphP[i].B[2], 2 ) );
    plot_slice();
    free( plot_info.data );
}

void gas_rho_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    plot_info.h = All.SofteningTable[0];
    plot_info.log_flag = 1;
    sprintf( plot_info.data_name, "rho" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( plot_info.ylabel, "" );
    sprintf( plot_info.title, "" );
    sprintf( plot_info.cb_label, "(10^x) gcm^{-3}Kpc^{-2}");
    plot_info.data = malloc( sizeof(double) * num );
    memset( plot_info.data, 0, sizeof(double) * num );
    plot_info.istart = SliceStart[0];
    plot_info.iend = SliceEnd[0];
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        plot_info.data[i-SliceStart[0]] = SphP[i].Density / ( cgs_g / CUBE(cgs_cm) );
    plot_slice();
    free( plot_info.data );
}

void mach_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    plot_info.h = All.SofteningTable[0];
    plot_info.log_flag = 1;
    sprintf( plot_info.data_name, "mn" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( plot_info.ylabel, "" );
    sprintf( plot_info.title, "" );
    sprintf( plot_info.cb_label, "(10^x)");
    plot_info.data = malloc( sizeof(double) * num );
    memset( plot_info.data, 0, sizeof(double) * num );
    plot_info.istart = SliceStart[0];
    plot_info.iend = SliceEnd[0];
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        plot_info.data[i-SliceStart[0]] = SphP[i].MachNumber;
    plot_slice();
    free( plot_info.data );
}

void divB_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    plot_info.h = All.SofteningTable[0];
    plot_info.log_flag = 1;
    sprintf( plot_info.data_name, "divB" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( plot_info.ylabel, "" );
    sprintf( plot_info.title, "" );
    sprintf( plot_info.cb_label, "(10^x)");
    plot_info.data = malloc( sizeof(double) * num );
    memset( plot_info.data, 0, sizeof(double) * num );
    plot_info.istart = SliceStart[0];
    plot_info.iend = SliceEnd[0];
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        plot_info.data[i-SliceStart[0]] = SphP[i].divB;
    plot_slice();
    free( plot_info.data );
}

void gas_vel_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    plot_info.h = All.SofteningTable[0];
    plot_info.log_flag = 1;
    sprintf( plot_info.data_name, "vel" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( plot_info.ylabel, "" );
    sprintf( plot_info.title, "" );
    sprintf( plot_info.cb_label, "(10^x)");
    plot_info.data = malloc( sizeof(double) * num );
    memset( plot_info.data, 0, sizeof(double) * num );
    plot_info.istart = SliceStart[0];
    plot_info.iend = SliceEnd[0];
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        plot_info.data[i-SliceStart[0]] = sqrt( pow( P[i].Vel[0], 2 ) +
                                                pow( P[i].Vel[1], 2 ) +
                                                pow( P[i].Vel[2], 2 ) );
    plot_slice();
    free( plot_info.data );
}

void hge_n_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    plot_info.h = All.SofteningTable[0];
    plot_info.log_flag = 1;
    plot_info.global_colorbar_flag = 1;
    sprintf( plot_info.data_name, "hge_n" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( plot_info.ylabel, "" );
    sprintf( plot_info.title, "" );
    sprintf( plot_info.cb_label, "(10^x)");
    plot_info.data = malloc( sizeof(double) * num );
    memset( plot_info.data, 0, sizeof(double) * num );
    plot_info.istart = SliceStart[0];
    plot_info.iend = SliceEnd[0];
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        plot_info.data[i-SliceStart[0]] = SphP[i].CRE_n0 * SphP[i].Density / cgs_mp * CUBE(cgs_cm);
    plot_slice();
    free( plot_info.data );
}

void cr_n_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    plot_info.h = All.SofteningTable[0];
    plot_info.log_flag = 1;
    plot_info.global_colorbar_flag = 1;
    sprintf( plot_info.data_name, "cr_n" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( plot_info.ylabel, "" );
    sprintf( plot_info.title, "" );
    sprintf( plot_info.cb_label, "(10^x)");
    plot_info.data = malloc( sizeof(double) * num );
    memset( plot_info.data, 0, sizeof(double) * num );
    plot_info.istart = SliceStart[0];
    plot_info.iend = SliceEnd[0];
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        plot_info.data[i-SliceStart[0]] = SphP[i].CR_n0 * SphP[i].Density  / cgs_me * CUBE(cgs_cm);
    plot_slice();
    free( plot_info.data );
}

void cr_e_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    plot_info.h = All.SofteningTable[0];
    plot_info.log_flag = 1;
    plot_info.global_colorbar_flag = 0;
    sprintf( plot_info.data_name, "cr_e" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( plot_info.ylabel, "" );
    sprintf( plot_info.title, "" );
    sprintf( plot_info.cb_label, "(10^x)");
    plot_info.data = malloc( sizeof(double) * num );
    memset( plot_info.data, 0, sizeof(double) * num );
    plot_info.istart = SliceStart[0];
    plot_info.iend = SliceEnd[0];
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        plot_info.data[i-SliceStart[0]] = SphP[i].CR_E0 * SphP[i].Density / cgs_erg * CUBE(cgs_cm);
    plot_slice();
    free( plot_info.data );
}

void hge_e_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    plot_info.h = All.SofteningTable[0];
    plot_info.log_flag = 1;
    plot_info.global_colorbar_flag = 0;
    sprintf( plot_info.data_name, "hge_e" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / All.MpcFlag );
    sprintf( plot_info.ylabel, "" );
    sprintf( plot_info.title, "" );
    sprintf( plot_info.cb_label, "(10^x)");
    plot_info.data = malloc( sizeof(double) * num );
    memset( plot_info.data, 0, sizeof(double) * num );
    plot_info.istart = SliceStart[0];
    plot_info.iend = SliceEnd[0];
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        plot_info.data[i-SliceStart[0]] = SphP[i].CRE_E0 * SphP[i].Density / cgs_erg * CUBE(cgs_cm);
    plot_slice();
    free( plot_info.data );
}

int compare_gas_rho( const void *a, const void *b ){
    return ((((struct sph_particle_data* )a)->Density) < (((struct sph_particle_data*)b)->Density)) ? 1: -1;
}

void sort_gas_rho(){
    int i;
    qsort( (void*)SphP, N_Gas, sizeof( struct sph_particle_data ), &compare_gas_rho );
    for ( i=0; i<10; i++ ) {
        printf( "%g\n", SphP[i].Density / ( cgs_g / CUBE(cgs_cm) ) );
    }
    printf( "\n" );
    for ( i=0; i<10; i++ ) {
        printf( "%g\n", SphP[N_Gas-10+i].Density / ( cgs_g / CUBE(cgs_cm) ) );
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

void analysis(){
    init_analysis();
    tree_build( 1 );
    tree_free();
    //magnetic_field_slice();
    //gas_rho_slice();
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
    free_analysis();
}
