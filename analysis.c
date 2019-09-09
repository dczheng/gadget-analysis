#include "allvars.h"

void init_analysis() {

    writelog( "initialize analysis...\n" );
    proj_k = All.ProjectDirection;
    proj_i = ( All.ProjectDirection + 1 ) % 3;
    proj_j = ( All.ProjectDirection + 2 ) % 3;
    Sproj = All.ProjectDirection + 'x';
    PicSize2 = SQR( All.PicSize );

    do_sync( "" );

    slice_init();
    if ( All.KernelInterpolation )
        init_kernel_matrix();
    init_img();

    create_dir( OutputDir );

    writelog( "initialize analysis... done.\n" );
    put_sep;

}

void free_analysis() {
    writelog( "free analysis ...\n" );
    if ( All.KernelInterpolation )
        free_kernel_matrix();
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

    if ( All.RadSpec ) {
        compute_particle_radio();
    }

    if ( ThisTask_Local == 0 ) {
        if ( All.TemperatureSlice ||
              All.Phase ||
              All.Group ||
              All.DensitySlice ||
              All.FieldCrenTDens ||
              All.HsmlTPdf ||
              All.UTPdf ||
              All.GasRatio ) {
            compute_temperature();
        }

        if ( All.Tree )
            tree_build();

        if ( All.FoF )
            fof();

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

        if ( All.CREeSlice ) {
            cree_slice();
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

        if ( All.MF )
            mass_function();

        if ( All.CorrGas ) {
            corr_gas();
        }
    
        if ( All.CorrDM ) {
            corr_dm();
        }
    
        if ( All.DensPdf ) {
            dens_pdf();
        }
        
        if ( All.CrenTPdf ) {
            cren_T_pdf();
        }
    
        if ( All.TPdf ) {
            T_pdf();
        }

        if ( All.GasRatio ) {
            gas_ratio();
        }

        if ( All.HsmlTPdf ) {
            hsml_T_pdf();
        }

        if ( All.HsmlDensPdf ) {
            hsml_dens_pdf();
        }

        if ( All.UTPdf ) {
            u_T_pdf();
        }

    }

    do_sync( "" );

    /*
    if ( All.RadSpec ) {
        compute_particle_radio();
    }
    */

    if ( ThisTask_Local == 0 && All.Group ) {
        sprintf( GroupDir, "%sgroup/", OutputDir );
        create_dir( GroupDir );
        group_analysis();
    }

    if ( All.TotSpec ) {
        total_radio_spectrum();
        do_sync( "total radio spectrum" );
    }

    if ( All.CrePressurePdf ) {
       cre_pressure_pdf();
    }

    if ( All.CorrTdiffDens ) {
        corr_Tdiff_dens();
    }

    if ( All.PdfTdiffDens ) {
        pdf_Tdiff_dens();
    }


    //tree_build();
    //tree_free();
    //output_rho();
    //vel_value();
    //

    if ( ThisTask_Local == 0 ) {

        if ( All.FoF )
            fof_free();
        if ( All.Tree )
            tree_free();

    }


    writelog( "analyais ... done.\n" );
    put_sep;


    free_analysis();
}
