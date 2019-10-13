#include "allvars.h"

void init_analysis() {

    int i;
    double u;

    put_header( "initialize analysis" );
    proj_k = All.ProjectDirection;
    proj_i = ( All.ProjectDirection + 1 ) % 3;
    proj_j = ( All.ProjectDirection + 2 ) % 3;
    Sproj = All.ProjectDirection + 'x';
    PicSize2 = SQR( All.PicSize );

    do_sync( "" );

    if ( All.KernelInterpolation )
        init_kernel_matrix();

    mymalloc1( ShortRangeTablePotential, sizeof(double) * NSRPTAB );

    for( i=0; i<NSRPTAB; i++ ) {
        u = 3.0 / NSRPTAB * ( i+0.5 );
        ShortRangeTablePotential[i] = erfc(u);
    }

    init_img();
    create_dir( OutputDir );
    put_end();

}

void free_analysis() {
    put_header( "free analysis" );
    if ( All.KernelInterpolation )
        free_kernel_matrix();
    free_img();

    myfree( ShortRangeTablePotential );

#ifdef TREE
    tree_free();
#endif
#ifdef FOF
    fof_free();
#endif


#ifndef ALTRAD
#ifdef RADSPEC
        free_particle_radio();
#endif
#endif
}

void analysis(){

    mytimer_start();
    put_header( "start analyais" );
    put_sep0;

    init_analysis();
    put_sep0;
    part_info();
    put_sep0;

    test_fof();
    test_group_pot();
#ifdef TREE
    tree_build();
    put_sep0;
#endif
#ifdef FOF
    fof();
    check_fof( 0, 1 );
    put_sep0;
#endif



#ifdef COMPUTETEMP
    compute_temperature();
#endif

#ifdef SMOOTH
    smooth();
    put_sep0;
#endif


#ifdef MACHNOISE
    remove_mach_noise();
#endif

#ifndef ALTRAD
#ifdef RADSPEC
    compute_particle_radio();
    put_sep0;
#endif
#endif

#ifdef TOTSPEC
    total_radio_spectrum();
    put_sep0;
#endif

#ifdef CREPPDF
    cre_pressure_pdf();
#endif

#ifdef CORRTDIFFDENS
    corr_Tdiff_dens();
#endif

#ifdef PDFTDIFFDENS
    pdf_Tdiff_dens();
#endif

    do_sync( "global compute" );
    put_sep0;

    if ( ThisTask_Local == 0 ) {
            slice();
#ifdef PHASE
            phase();
#endif

#ifdef BPDF
            B_Pdf();
#endif

#ifdef DIVBERRPDF
            DivB_Err_Pdf();
#endif

#ifdef POWSPEC
            powerspec();
#endif

#ifdef DENSITYSLICE
            density_slice();
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

#ifdef MACHDENSPDF
            mach_dens_pdf();
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

#ifdef DIVBERRDENSPDF
            divBerr_dens_pdf();
#endif

    }

    do_sync( "" );
    put_sep0;

    free_analysis();
    mytimer_end();
    put_end();
}
