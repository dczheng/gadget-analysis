#include "allvars.h"

void init_analysis() {

    writelog( "initialize analysis...\n" );
    All.proj_k = All.ProjectDirection;
    All.proj_i = ( All.ProjectDirection + 1 ) % 3;
    All.proj_j = ( All.ProjectDirection + 2 ) % 3;
    All.Sproj = All.ProjectDirection + 'x';
    All.PicSize2 = SQR( All.PicSize );

    if ( All.RadSpec ) {
        compute_particle_radio();
    }
    do_sync( "" );

    slice_init();
    if ( All.KernelInterpolation )
        init_kernel_matrix();
    init_conv_kernel();
    init_img();

    create_dir( All.OutputDir );

    writelog( "initialize analysis... done.\n" );
    put_sep;

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
    put_sep;
}

void analysis(){

    init_analysis();

    writelog( "analyais ...\n" );

    put_sep0;
    //
    part_info();

//    endrun(20190625);

    if ( ThisTask_Local == 0 ) {
        if ( (All.TemperatureSlice ||
              All.Phase ||
              All.Group )
             && All.ReadTemp == 0  ) {
            compute_temperature();
        }

        if ( All.Phase ) {
            phase();
        }

        if ( All.MachSlice ) {
            mach_slice();
        }

        if ( All.BSlice ) {
            mag_slice();
        }

        if ( All.CREnSlice ) {
            cren_slice();
        }

        if ( All.BPdf ) {
            B_Pdf();
        }

        if ( All.RadSlice ) {
            radio_slice();
        }

        if ( All.PowSpec )
            powerspec();

        if ( All.DensitySlice )
            density_slice();

        if ( All.TemperatureSlice )
            temperature_slice();

        if ( All.FoF )
            fof();

        if ( All.MF )
            mass_function();

    }

    do_sync( "" );

    /*
    if ( All.RadSpec ) {
        compute_particle_radio();
    }
    */

    if ( ThisTask_Local == 0 && All.Group ) {
        sprintf( All.GroupDir, "%sgroup/", All.OutputDir );
        create_dir( All.GroupDir );
        group_analysis();
    }

    if ( All.TotSpec ) {
        total_radio_spectrum();
        do_sync( "total radio spectrum" );
    }

    if ( ThisTask_Local == 0 && All.FoF ) {
        fof_free();
    }

    if ( All.CrePressurePdf ) {
       cre_pressure_pdf();
    }

    //tree_build();
    //tree_free();
    //output_rho();
    //vel_value();
    //
    writelog( "analyais ... done.\n" );
    put_sep;

    free_analysis();
}
