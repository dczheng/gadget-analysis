#include "allvars.h"

#define F( c, a, q1, q2, p )                 ( ( p<q1 || p>q2 ) ? 0 : c * pow( p, -a ) )

#define BETA(a, q)                           beta( (a-2)*0.5, (3-a)*0.5, 1/(1+q*q) )

#define NUMBER_DENSITY( c, a, q )            ( c * pow( q, 1.0-a ) / ( a-1.0 ) )

#define ENERGY( c, a, q )                    ( guc.c2 * c / (a-1.0) * ( 0.5 * BETA(a, q) + (sqrt( 1+q*q ) - 1.0 ) * pow( q, 1-a ) ) )
#define ENERGY2( c, a, q1, q2 )              ( ENERGY0( c, a, q1 ) - ENERGY0( c, a, q2 ) )

#define MEAN_ENERGY( a, q )                  ( ( 0.5 * pow( q, a-1.0 ) * BETA(a, q) + sqrt( 1+q*q ) - 1.0 ) * guc.mec2 )
#define MEAN_ENERGY2( a, q1, q2 )            ( ENERGY02( 1, a, q1, q2 ) / NUMBER_DENSITY2( 1, a, q1, q2 ) * guc.m_e )

#define PRESSURE( c, a, q )                  ( guc.c2 * c / 6 * BETA( a, q ) )
#define PRESSURE2( c, a, q1, q2 )            ( PRESSURE( c, a, q1 ) - PRESSURE( c, a, q2 ) )

#define hge_pressure( i )                    ( PRESSURE2( SphP[i].CRE_C, SphP[i].CRE_Alpha, SphP[i].CRE_qmin, SphP[i].CRE_qmax ) )


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

void compute_hge_pressure() {

    long i, t;
    t = N_Gas / 10;
    put_sep0;
    mytimer_start();
    writelog( "compute hge pressure ...\n" );

    for( i=0; i<N_Gas; i++ ) {
        if ( i % t == 0 )
            writelog( "[%10li] [%10li] [%5.1f%%]\n", i, N_Gas, ( (double)(i) ) / N_Gas * 100 );
        if ( i % NTask_Local != ThisTask_Local )
            continue;
        if ( SphP[i].CRE_n != 0 ) {
            SphP[i].CRE_P = hge_pressure( i ) * SphP[i].Density;
        }
    }
    mytimer_end();
    put_sep0;

}

void hge_pressure_pdf() {

    long i;
    int PicSize, ii, jj;
    char buf[100];
    double p, hge_r, cr_r, hge_r_max, cr_r_max, hge_r_min, cr_r_min,
           dlogcr_r, dloghge_r, *img, sum;

    PicSize = All.PicSize;
    hge_r_max = cr_r_max = 0;
    hge_r_min = cr_r_min = DBL_MAX;

    compute_hge_pressure();
    do_sync_local( "" );

    if ( ThisTask_Local != 0 )
        return;

    mytimer_start();
    writelog( "hge pressure ...\n" );

    for( i=0; i<N_Gas; i++ ) {

         if ( SphP[i].CRE_n == 0 )
             continue;

        p = get_pressure( i );

        hge_r = SphP[i].CRE_P / p;
        cr_r  = SphP[i].CR_P0 / p;

        vmax2( hge_r_max, hge_r );
        vmax2( cr_r_max, cr_r );

        vmin2( hge_r_min, hge_r, 0);
        vmin2( cr_r_min, cr_r, 0 );

    }

    writelog( "[min] cr: %g, hge: %g\n", cr_r_min, hge_r_min );
    writelog( "[max] cr: %g, hge: %g\n", cr_r_max, hge_r_max );

    dlogcr_r = log10( cr_r_max / cr_r_min ) / ( PicSize-1 );
    dloghge_r = log10( cr_r_max / cr_r_min ) / ( PicSize-1 );

    mymalloc2( img, sizeof(double) * SQR(PicSize) );

    sum = 0;
    for( i=0; i<N_Gas; i++ ) {

        if ( SphP[i].CRE_n == 0 )
            continue;

        p = get_pressure( i );
        hge_r = SphP[i].CRE_P / p;
        cr_r  = SphP[i].CR_P0 / p;

        ii = ( log10( hge_r / hge_r_min ) ) / dloghge_r;
        jj = ( log10( cr_r / cr_r_min ) ) / dlogcr_r;

        check_picture_index( ii );
        check_picture_index( jj );

        img[ ii*PicSize + jj ] += 1;
        sum ++;

    }

    for( i=0; i<SQR(PicSize); i++ )
        img[i] /= sum * dloghge_r * dlogcr_r;

    create_dir( "./HgePressurePdf" );
    sprintf( buf, "./HgePressurePdf/HgePressurePdf_%.2f.dat", All.RedShift );

    image.img = img;
    img_xmin = log10( cr_r_min );
    img_xmax = log10( cr_r_max );
    img_ymin = log10( hge_r_min );
    img_ymax = log10( hge_r_max );

    write_img1( buf, "HgePressurePdf" );

    myfree( img );

    mytimer_end();
    put_sep0;

}
