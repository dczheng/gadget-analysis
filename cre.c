#include "allvars.h"

#define F( c, a, q1, q2, p )                 ( ( p<q1 || p>q2 ) ? 0 : c * pow( p, -a ) )

#define BETA(a, q)                           beta( (a-2)*0.5, (3-a)*0.5, 1/(1+q*q) )

#define NUMBER_DENSITY( c, a, q )            ( (c) * pow( (q), 1.0-(a) ) / ( (a)-1.0 ) )
#define NUMBER_DENSITY2( c, a, q1, q2 )      ( NUMBER_DENSITY( (c), (a), (q1) ) - NUMBER_DENSITY( (c), (a), (q2) ) )

#define ENERGY( c, a, q )                    ( guc.c2 * (c) / ((a)-1.0) * ( 0.5 *  BETA((a), (q)) + (sqrt( 1+(q)*(q)) - 1.0 ) * pow( (q), 1-(a) ) ) )
#define ENERGY2( c, a, q1, q2 )              ( ENERGY( (c), (a), (q1) ) -  ENERGY( (c), (a), (q2) ) )

#define MEAN_ENERGY( a, q )                  ( ( 0.5 * pow( (q), (a)-1.0 ) * BETA((a), (q)) + sqrt( 1+((q)*(q)) ) - 1.0 ) * guc.mec2 )
#define MEAN_ENERGY2( a, q1, q2 )            ( ENERGY2( 1, (a), (q1), (q2) ) / NUMBER_DENSITY2( 1, (a), (q1), (q2) ) * guc.m_e )


#define PRESSURE( c, a, q )                  ( guc.c2 * (c) / 6 * BETA( (a), (q) ) )
#define PRESSURE2( c, a, q1, q2 )            ( PRESSURE( (c), (a), (q1) ) - PRESSURE( (c), (a), (q2) ) )

#define cre_pressure( i )                    ( PRESSURE2( SphP[i].CRE_C, SphP[i].CRE_Alpha, SphP[i].CRE_qmin, SphP[i].CRE_qmax ) )


double beta_inte( double x, void *params ) {

    double *p = params;
    return pow( x, p[0]-1.0 ) * pow( 1.0-x, p[1]-1.0 );

}

double beta( double a, double b, double x ) {

    double r, err;
    gsl_function F;
    double params[2];

    /*
    gsl_integration_workspace *inte_ws_local;
    inte_ws_local = gsl_integration_workspace_alloc( GSL_INTE_WS_LEN );
    */

    if ( a<0 ) {
        RAISE_SIGSTOP();
        endrun( 20190922 );
    }

    F.function = &beta_inte;

    if ( b>0.0 )
        return gsl_sf_beta_inc( a, b, x ) * gsl_sf_beta( a, b );

    F.params = params;
    params[0] = a;
    params[1] = b;

    gsl_integration_qag( &F, 0, x,
            GSL_INTE_ERR_ABS, GSL_INTE_ERR_REL, GSL_INTE_WS_LEN, GSL_INTE_KEY,
            inte_ws, &r, &err );

    //gsl_integration_workspace_free( inte_ws_local );

    return r;

}

void compute_cre_pressure() {

    long i, t;
    t = N_Gas / 10;
    mytimer_start();
    writelog( "compute cre pressure ...\n" );

    for( i=0; i<N_Gas; i++ ) {
        if ( i % t == 0 )
            writelog( "[%10li] [%10li] [%5.1f%%]\n", i, N_Gas, ( (double)(i) ) / N_Gas * 100 );
        if ( i % NTask_Local != ThisTask_Local )
            continue;
        if ( SphP[i].CRE_n != 0 ) {
            SphP[i].CRE_P = cre_pressure( i ) * SphP[i].Density;
        }
    }
    mytimer_end();

}

void cre_pressure_pdf() {

    long i;
    int PicSize, ii, jj, *cre_cr, *cr_bar;
    char buf[100], bins=30;
    double p, cre_r, cr_r, cre_r_max, cr_r_max, cre_r_min, cr_r_min,
           dlogcr_r, dlogcre_r, sum,
           dlogcre_bar, dlogcr_bar, t;

    FILE *fd;

    PicSize = All.PicSize;
    cre_r_max = cr_r_max = 0;
    cre_r_min = cr_r_min = DBL_MAX;

    compute_cre_pressure();
    do_sync_local( "" );

    if ( ThisTask_Local != 0 )
        return;

    mytimer_start();
    writelog( "cre pressure ...\n" );

    for( i=0; i<N_Gas; i++ ) {

         if ( SphP[i].CRE_n == 0 )
             continue;

        p = get_pressure( i );

        cre_r = SphP[i].CRE_P / p;
        cr_r  = SphP[i].CR_P0 / p;

        //printf( "CRE_P: %g, CR_P: %g\n", SphP[i].CRE_P, SphP[i].CR_P0 );

        vmax2( cre_r_max, cre_r );
        vmax2( cr_r_max, cr_r );

        vmin2( cre_r_min, cre_r );
        vmin2( cr_r_min, cr_r );

    }

    writelog( "[min] cr: %g, cre: %g\n", cr_r_min, cre_r_min );
    writelog( "[max] cr: %g, cre: %g\n", cr_r_max, cre_r_max );

    dlogcr_r = log10( cr_r_max / cr_r_min ) / ( PicSize-1 );
    dlogcre_r = log10( cre_r_max / cre_r_min ) / ( PicSize-1 );

    reset_img();

    sum = 0;
    for( i=0; i<N_Gas; i++ ) {

        if ( SphP[i].CRE_n == 0 )
            continue;

        p = get_pressure( i );
        cre_r = SphP[i].CRE_P / p;
        cr_r  = SphP[i].CR_P0 / p;

        ii = ( log10( cre_r / cre_r_min ) ) / dlogcre_r;
        jj = ( log10( cr_r / cr_r_min ) ) / dlogcr_r;

        check_picture_index( ii );
        check_picture_index( jj );

        image.img[ ii*PicSize + jj ] += 1;
        sum ++;

    }


    for( i=0; i<SQR(PicSize); i++ )
        image.img[i] /= sum * dlogcre_r * dlogcr_r;

    create_dir( "%s/CrePressurePdf", OutputDir);
    sprintf( buf, "%s/CrePressurePdf/CrePressurePdf_%03i.dat",
                    OutputDir,  SnapIndex );

    img_xmin = log10( cr_r_min );
    img_xmax = log10( cr_r_max );
    img_ymin = log10( cre_r_min );
    img_ymax = log10( cre_r_max );

    write_img( buf, "CrePressurePdf" );

    dlogcre_bar = log10( cre_r_max / cre_r_min ) / ( bins-1 );
    dlogcr_bar = log10( cr_r_max / cr_r_min ) / ( bins-1 );

    mymalloc2( cre_cr, sizeof(int) * bins );
    mymalloc2( cr_bar, sizeof(int) * bins );

    for( i=0; i<N_Gas; i++ ){

        if ( SphP[i].CRE_n == 0 )
            continue;

        p = get_pressure( i );
        t = SphP[i].CRE_P / p;
        ii = log10(t/cre_r_min) / dlogcre_bar;
        if ( ii >= bins )
            printf( "%i\n", ii );
        cre_cr[ii] ++;

        t = SphP[i].CR_P0 / p;
        ii = log10(t/cr_r_min) / dlogcr_bar;
        if ( ii >= bins )
            printf( "%i\n", ii );
        cr_bar[ii] ++;
    }

    fd = myfopen( "w", "%s/CrePressurePdf/CrePressure_%03i.dat",
                        OutputDir, SnapIndex );

    for ( i=0; i<bins; i++ ){
        fprintf( fd, "%g %i %g %i\n",
                pow( 10, log10(cre_r_min) + dlogcre_bar * i ),
                cre_cr[i],
                pow( 10, log10(cr_r_min) + dlogcr_bar * i ),
                cr_bar[i]
                );
    }

    fclose( fd );
    myfree( cre_cr );
    myfree( cr_bar );
    mytimer_end();

}

