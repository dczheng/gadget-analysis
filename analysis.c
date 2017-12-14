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
    plot_info.h = para.SofteningTable[0];
    plot_info.log_flag = 1;
    sprintf( plot_info.data_name, "B" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / para.MpcFlag );
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

    plot_info.h = para.SofteningTable[0];
    plot_info.log_flag = 1;
    sprintf( plot_info.data_name, "rho" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / para.MpcFlag );
    sprintf( plot_info.ylabel, "" );
    sprintf( plot_info.title, "" );
    sprintf( plot_info.cb_label, "(10^x)");
    plot_info.data = malloc( sizeof(double) * num );
    memset( plot_info.data, 0, sizeof(double) * num );
    plot_info.istart = SliceStart[0];
    plot_info.iend = SliceEnd[0];
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ )
        plot_info.data[i-SliceStart[0]] = SphP[i].Density / ( cgs_g/(CUBE(cgs_cm)) );
    plot_slice();
    free( plot_info.data );
}

void mach_slice() {
    int num, i;
    num = SliceEnd[0] - SliceStart[0];

    plot_info.h = para.SofteningTable[0];
    plot_info.log_flag = 1;
    sprintf( plot_info.data_name, "mn" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / para.MpcFlag );
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

    plot_info.h = para.SofteningTable[0];
    plot_info.log_flag = 1;
    sprintf( plot_info.data_name, "divB" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / para.MpcFlag );
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

    plot_info.h = para.SofteningTable[0];
    plot_info.log_flag = 1;
    sprintf( plot_info.data_name, "vel" );
    sprintf( plot_info.xlabel, "%g Mpc", proj_size / para.MpcFlag );
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

void analysis(){
    init_analysis();
    magnetic_field_slice();
    gas_rho_slice();
    divB_slice();
    gas_vel_slice();
    mach_slice();
    free_analysis();
}
