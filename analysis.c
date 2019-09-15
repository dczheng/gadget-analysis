#include "allvars.h"

void init_analysis() {

    put_header( "initialize analysis" );
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
}

void free_analysis() {
    put_header( "free analysis" );
    if ( All.KernelInterpolation )
        free_kernel_matrix();
    free_img();

#ifdef RADSPEC
        free_particle_radio();
#endif
}

void analysis(){


    mytimer_start();
    put_header( "start analyais" );
    put_sep0;
    init_analysis();
    put_sep0;

    if ( ThisTask_Local == 0 ) {

#ifdef TREE
            tree_build();
            put_sep0;
#endif

#ifdef SMOOTH
        smooth();
        put_sep0;
#endif

#ifdef FOF
            fof();
            put_sep0;
#endif

#ifdef COMPUTETEMP
            compute_temperature();
#endif

    }

    do_sync( "global compute1" );

#ifdef RADSPEC
        compute_particle_radio();
        put_sep0;
#endif

#ifdef TOTSPEC
        total_radio_spectrum();
        put_sep0;
#endif

    do_sync( "global compute2" );
    put_sep0;

    part_info();
    put_sep0;

#ifdef CREPPDF
       cre_pressure_pdf();
#endif

#ifdef CORRTDIFFDENS
        corr_Tdiff_dens();
#endif

#ifdef PDFTDIFFDENS
        pdf_Tdiff_dens();
#endif

    if ( ThisTask_Local == 0 ) {
#ifdef PHASE
            phase();
#endif

#ifdef MACHSLICE
            mach_slice();
#endif

#ifdef BSLICE
            mag_slice();
#endif

#ifdef CRENSLICE
            cren_slice();
#endif

#ifdef CREESLICE
            cree_slice();
#endif

#ifdef BPDF
            B_Pdf();
#endif

#ifdef DIVBERRPDF
            DivB_Err_Pdf();
#endif

#ifdef RADSLICE
            radio_slice();
#endif

#ifdef POWSPEC
            powerspec();
#endif

#ifdef DENSITYSLICE
            density_slice();
#endif

#ifdef TEMPSLICE
            temperature_slice();
#endif

#ifdef MF
            mass_function();
#endif

#ifdef CORRGAS
            corr_gas();
#endif
    
#ifdef CORRDM
            corr_dm();
#endif
    
#ifdef DENSPDF
            dens_pdf();
#endif
        
#ifdef CRENTPDF
            cren_T_pdf();
#endif
    
#ifdef TPDF
            T_pdf();
#endif

#ifdef GASRATIO
            gas_ratio();
#endif

#ifdef HSMLTPDF
            hsml_T_pdf();
#endif

#ifdef HSMLDENSPDF
            hsml_dens_pdf();
#endif

#ifdef UTPDF
            u_T_pdf();
#endif

#ifdef GROUP
            sprintf( GroupDir, "%sgroup/", OutputDir );
            create_dir( GroupDir );
            group_analysis();
#endif

#ifdef BDENSPDF
            B_dens_pdf();
#endif

    }
    do_sync( "" );
    put_sep0;

    if ( ThisTask_Local == 0 ) {

#ifdef FOF
            fof_free();
#endif
#ifdef TREE
            tree_free();
#endif
    }

    free_analysis();
    mytimer_end();
}
