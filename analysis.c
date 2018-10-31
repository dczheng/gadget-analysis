#include "allvars.h"

void init_analysis() {
    writelog( "initialize analysis...\n" );
    All.proj_k = All.ProjectDirection;
    All.proj_i = ( All.ProjectDirection + 1 ) % 3;
    All.proj_j = ( All.ProjectDirection + 2 ) % 3;
    All.Sproj = All.ProjectDirection + 'x';
    All.PicSize2 = SQR( All.PicSize );
    slice();
    if ( All.KernelInterpolation )
        init_kernel_matrix();
    init_conv_kernel();
    init_img();

    writelog( "initialize analysis... done.\n" );
    put_block_line;
}

void free_analysis() {
    writelog( "free analysis ...\n" );
    if ( All.KernelInterpolation )
        free_kernel_matrix();
    free_conv_kernel();
    free_img();

    if ( All.RadSpec ) {
        free_particle_radio();
    }

    writelog( "free analysis ... done.\n" );
    put_block_line;
}


void analysis(){

    init_analysis();
    //printf( "%g\n", All.RedShift );
    if ( (All.GasTemperature ||
          All.GasState ||
          All.Group )
         && All.ReadTemp == 0  ) {
        compute_temperature();
    }

    if ( All.GasState ) {
        gas_state();
    }

    if ( All.PowSpec )
        powerspec();

    if ( All.GasDensity )
        gas_density();

    if ( All.GasTemperature )
        gas_temperature();

    if ( All.FoF ) {
        fof();
        do_sync( "fof" );
    }


    if ( All.RadSpec )
        compute_particle_radio();

    if ( All.MF )
        mass_function();


    if ( All.Group ) {
        create_dir( All.GroupDir );
        group_analysis();

    }

    if ( All.TotSpec ) {
        total_radio_spectrum();
        do_sync( "total radio spectrum" );
    }

    if( All.FoF )
        fof_free();




    //tree_build();
    //tree_free();
    //output_rho();
    //vel_value();
    //sort_gas_rho();
    free_analysis();
}
