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
}

void analysis(){


    mytimer_start();
    writelog( "start analyais ...\n" );
    put_sep0;
    init_analysis();
    put_sep0;

    if ( ThisTask_Local == 0 ) {
        if ( All.Tree ) {
            tree_build();
            put_sep0;
        }

        smooth();

        if ( All.FoF ) {
            fof();
            put_sep0;
        }

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

    }

    do_sync("");

    if ( All.RadSpec ) {
        compute_particle_radio();
        put_sep0;
    }

    if ( All.TotSpec ) {
        total_radio_spectrum();
        put_sep0;
    }

    do_sync( "global compute" );
    put_sep0;

    part_info();
    put_sep0;

    if ( All.CrePressurePdf ) {
       cre_pressure_pdf();
    }

    if ( All.CorrTdiffDens ) {
        corr_Tdiff_dens();
    }

    if ( All.PdfTdiffDens ) {
        pdf_Tdiff_dens();
    }

    if ( ThisTask_Local == 0 ) {
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

        if ( All.PowSpec ) {
            powerspec();
        }

        if ( All.DensitySlice ) {
            density_slice();
        }

        if ( All.TemperatureSlice ) {
            temperature_slice();
        }

        if ( All.MF ) {
            mass_function();
        }

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

        if ( All.Group ) {
            sprintf( GroupDir, "%sgroup/", OutputDir );
            create_dir( GroupDir );
            group_analysis();
        }

    }
    do_sync( "" );
    put_sep0;

    if ( ThisTask_Local == 0 ) {

        if ( All.FoF ) {
            fof_free();
        }
        if ( All.Tree ) {
            tree_free();
        }
    }

    writelog( "analyais ... done.\n" );
    free_analysis();
    mytimer_end();
}
