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

    slice_init();
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

    tree_free();
    fof_free();
    free_particle_radio();
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

    tree_build();

    fof();

    compute_temperature();
    smooth();
    remove_mach_noise();
    compute_particle_radio();
    total_radio_spectrum();

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
            density_slice();
            phase();
            B_Pdf();
            DivB_Err_Pdf();
            powerspec();
            mass_function();

            corr_gas();
            corr_dm();
    
            dens_pdf();
            cren_T_pdf();
            T_pdf();
            hsml_T_pdf();
            gas_ratio();
            hsml_dens_pdf();
            mach_dens_pdf();
            u_T_pdf();
            B_dens_pdf();
            divBerr_dens_pdf();

    }

    sprintf( GroupDir, "%sgroup/", OutputDir );
    create_dir( GroupDir );
    group_analysis();


    do_sync( "" );
    put_sep0;

    free_analysis();
    mytimer_end();
    put_end();
}
