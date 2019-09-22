#include"allvars.h"


void pdf2d_or_field2d( double *x, double *y, double *w, long num, char *dn, 
             int flag, double *mm, int Nmin ) {
    /*
    flag:
        bit 1: mode, 0 for pdf, 1 for 2d field 
        bit 2: x log flag
        bit 3: y log flag
        bit 4: fixed x range 
        bit 5: fixed y range
        bit 6: nomalize
    require data > 0 for log.
    */

    double dx, dy, xmin, xmax, ymin, ymax, s, area;
    int i, j, N, N2;
    long p;

    char buf[100];
    reset_img();

    N = PicSize;
    N2 = SQR(N);

    if ( (flag & PDF2D_BIT_FIXEDX ) == 0  ||
         (flag & PDF2D_BIT_FIXEDY ) == 0 ) {
        
        xmin = ymin = DBL_MAX;
        xmax = ymax = -DBL_MAX;
        for( p=0; p<num; p++ ) {
            vmin2( xmin, x[p] );
            vmax2( xmax, x[p] );
            vmin2( ymin, y[p] );
            vmax2( ymax, y[p] );
        }
    }

    if ( flag & PDF2D_BIT_FIXEDX ) {
        xmin = mm[0];
        xmax = mm[1];
    }
    if ( flag & PDF2D_BIT_FIXEDY ) {
        ymin = mm[2];
        ymax = mm[3];
    }

    if ( flag & PDF2D_BIT_XLOG ) {
        if ( xmin <= 0 || xmax <= 0 )
            endruns( "log scale requires positive value");
        xmin = log10( xmin );
        xmax = log10( xmax );
        for( p=0; p<num; p++ )
            x[p] = log10(x[p]);
    }

    if ( flag & PDF2D_BIT_YLOG ) {
        if ( ymin <= 0 || ymax <= 0 )
            endruns( "log scale requires positive value");
        ymin = log10( ymin );
        ymax = log10( ymax );
        for( p=0; p<num; p++ )
            y[p] = log10(y[p]);
    }

    dx = ( xmax - xmin ) / ( N-1 );
    dy = ( ymax - ymin ) / ( N-1 );
    area = dx * dy;

    writelog( "[%s], dx: %g ,dy: %g, area: %g\n", __FUNCTION__,
        dx, dy, area );

    for( p=0; p<N_Gas; p++ ) {
        i = ( y[p] - ymin ) / dy;
        j = ( x[p] - xmin ) / dx;
        if ( i < 0 || i > N-1 || j < 0 || j > N-1 )
            continue;

        if ( flag & PDF2D_BIT_MODE  ) {
            if ( NULL == w )
                endruns( "require w data" );
            image.img[ i*N + j ] += w[p];
            image.num[ i*N + j ] ++;

        }
        else {
            if ( w )
                image.img[ i*N + j ] += w[p];
            else
                image.img[ i*N + j ] ++;
        }

    }

    if ( flag & PDF2D_BIT_MODE  ) {
        for( p=0; p<N2; p++ )
            if ( image.num[p] >= Nmin )
                image.img[p] /= image.num[p];
            else
                image.img[p] = 0;
    }
    else {
        if ( flag & PDF2D_BIT_UNIT_AREA )
            for( p=0; p<N2; p++ )
                image.img[p] /=  area;

        if ( flag & PDF2D_BIT_NORM ) {
            for( s=0,p=0; p<N2; p++ )
                s += image.img[p];
            for( p=0; p<N2; p++ )
                image.img[p] /= s;
        }
    }

    create_dir( "%s%s", OutputDir, dn );
    sprintf( buf, "%s%s/%s_%03i.dat", OutputDir, dn, dn, SnapIndex );

    img_xmin = xmin;
    img_xmax = xmax;
    img_ymin = ymin;
    img_ymax = ymax;

    if ( flag & PDF2D_BIT_XLOG )
        img_xlog = 1;
    else
        img_xlog = 0;

    if ( flag & PDF2D_BIT_YLOG )
        img_ylog = 1;
    else
        img_ylog = 0;

    write_img( buf, dn );

}

#ifdef BDENSPDF
void B_dens_pdf() {

    double *bdp_x, *bdp_y, mm[4];
    long p; 
    int flag;

    put_header( "B-Density-Pdf" );
    mymalloc1( bdp_x, sizeof(double) * N_Gas );
    mymalloc1( bdp_y, sizeof(double) * N_Gas );

    for( p=0; p<N_Gas; p++ ) {
        bdp_x[p] = SphP[p].Density / Time3 / RhoBaryon;
        bdp_y[p] = get_B( p ) * 1e6;
    }

    flag = 0;
    flag |= PDF2D_BIT_XLOG;
    flag |= PDF2D_BIT_YLOG;

    if ( All.BDensPdfDensMin > 0 || 
            All.BDensPdfDensMax > 0 ) {
        if ( All.BDensPdfDensMin * All.BDensPdfDensMax == 0 )
            endrun( 20180810 );

        flag |= PDF2D_BIT_FIXEDY;
        mm[0] = All.BDensPdfDensMin;
        mm[1] = All.BDensPdfDensMax;

    }

    if ( All.BDensPdfBMin > 0 || 
            All.BDensPdfBMax > 0 ) {
        if ( All.BDensPdfBMin * All.BDensPdfBMax == 0 )
            endrun( 20180810 );

        flag |= PDF2D_BIT_FIXEDX;
        mm[2] = All.BDensPdfBMin*1e6;
        mm[3] = All.BDensPdfBMax*1e6;

    }

    pdf2d( bdp_x, bdp_y, NULL, N_Gas, "BDensPdf", flag, mm );

    myfree( bdp_x );
    myfree( bdp_y );
    put_end();
}
#endif

#ifdef CRENTPDF
void cren_T_pdf() {

    double *ntp_x, *ntp_y, mm[4], n;
    long p, N; 
    int flag;

    mymalloc1( ntp_x, sizeof(double) * N_Gas );
    mymalloc1( ntp_y, sizeof(double) * N_Gas );
    for( p=0, N=0; p<N_Gas; p++ ) {
        n = SphP[p].CRE_n * SphP[p].Density / guc.m_e * ( 1/CUBE(g2c.cm) ) ;
        if ( n == 0 )
            continue;

        ntp_x[N] = n;
        ntp_y[N] = SphP[p].Temp;
        N++;
    }


    flag = 0;
    flag |= PDF2D_BIT_XLOG;
    flag |= PDF2D_BIT_YLOG;

    if ( All.CrenTPdfnMin > 0 || 
            All.CrenTPdfnMax > 0 ) {
        if ( All.CrenTPdfnMin * All.CrenTPdfnMax == 0 )
            endrun( 20180810 );

        flag |= PDF2D_BIT_FIXEDX;
        mm[0] = All.CrenTPdfnMin;
        mm[1] = All.CrenTPdfnMax;

    }

    if ( All.CrenTPdfTMin > 0 || 
            All.CrenTPdfTMax > 0 ) {
        if ( All.CrenTPdfTMin * All.CrenTPdfTMax == 0 )
            endrun( 20180810 );

        flag |= PDF2D_BIT_FIXEDY;
        mm[2] = All.CrenTPdfTMin;
        mm[3] = All.CrenTPdfTMax;

    }

    pdf2d( ntp_x, ntp_y, NULL, N, "CrenTPdf", flag, mm );

    myfree( ntp_x );
    myfree( ntp_y );
}
#endif

#ifdef HSMLTPDF
void hsml_T_pdf() {

    double *htp_x, *htp_y;
    long p; 
    int flag;

    mymalloc1( htp_x, sizeof(double) * N_Gas );
    mymalloc1( htp_y, sizeof(double) * N_Gas );
    for( p=0; p<N_Gas; p++ ) {
        htp_x[p] = SphP[p].Hsml;
        /*
        printf( "%g %g %g\n", SphP[p].Hsml, SphP[p].Density,
            SphP[p].Hsml*SphP[p].Density );
        */
        htp_y[p] = SphP[p].Temp;
    }

    flag = 0;
    flag |= PDF2D_BIT_XLOG;
    flag |= PDF2D_BIT_YLOG;

    pdf2d( htp_x, htp_y, NULL, N_Gas, "HsmlTPdf", flag, NULL );

    myfree( htp_x );
    myfree( htp_y );
}
#endif

#ifdef UTPDF
void u_T_pdf() {

    double *utp_x, *utp_y;
    long p; 
    int flag;

    mymalloc1( utp_x, sizeof(double) * N_Gas );
    mymalloc1( utp_y, sizeof(double) * N_Gas );
    for( p=0; p<N_Gas; p++ ) {
        utp_x[p] = SphP[p].u;
        utp_y[p] = SphP[p].Temp;
    }

    flag = 0;
    flag |= PDF2D_BIT_XLOG;
    flag |= PDF2D_BIT_YLOG;

    pdf2d( utp_x, utp_y, NULL, N_Gas, "UTPdf", flag, NULL );

    myfree( utp_x );
    myfree( utp_y );
}
#endif

#ifdef HSMLDENSPDF
void hsml_dens_pdf() {

    double *htp_x, *htp_y;
    long p; 
    int flag;

    mymalloc1( htp_x, sizeof(double) * N_Gas );
    mymalloc1( htp_y, sizeof(double) * N_Gas );
    for( p=0; p<N_Gas; p++ ) {
        htp_x[p] = SphP[p].Hsml;
        /*
        printf( "%g %g %g\n", SphP[p].Hsml, SphP[p].Density,
            SphP[p].Hsml*SphP[p].Density );
        */
        htp_y[p] = SphP[p].Density / Time3 / RhoBaryon;
    }

    flag = 0;
    flag |= PDF2D_BIT_XLOG;
    flag |= PDF2D_BIT_YLOG;

    pdf2d( htp_x, htp_y, NULL, N_Gas, "HsmlDensPdf", flag, NULL );

    myfree( htp_x );
    myfree( htp_y );
}
#endif

#ifdef DIVBERRDENSPDF
void divBerr_dens_pdf() {

    double *divbdp_x, *divbdp_y, mm[4], B;
    long p, N; 
    int flag;

    put_header( "divBerr-Density-Pdf" );
    mymalloc1( divbdp_x, sizeof(double) * N_Gas );
    mymalloc1( divbdp_y, sizeof(double) * N_Gas );

    for( p=0, N=0; p<N_Gas; p++ ) {
        B = get_B(p);
        if ( B<=0 || SphP[p].divB == 0 )
            continue;
        divbdp_x[N] = SphP[p].Density / Time3 / RhoBaryon;
        divbdp_y[N] = fabs(SphP[p].divB) * SphP[p].Hsml / B;
        N++;
    }

    flag = 0;
    flag |= PDF2D_BIT_XLOG;
    flag |= PDF2D_BIT_YLOG;

    if ( All.DivBerrDensPdfDensMin > 0 || 
            All.DivBerrDensPdfDensMax > 0 ) {
        if ( All.DivBerrDensPdfDensMin * All.DivBerrDensPdfDensMax == 0 )
            endrun( 20180810 );

        flag |= PDF2D_BIT_FIXEDY;
        mm[0] = All.DivBerrDensPdfDensMin;
        mm[1] = All.DivBerrDensPdfDensMax;

    }

    if ( All.DivBerrDensPdfDivBMin > 0 || 
            All.DivBerrDensPdfDivBMax > 0 ) {
        if ( All.DivBerrDensPdfDivBMin * All.DivBerrDensPdfDivBMax == 0 )
            endrun( 20180810 );

        flag |= PDF2D_BIT_FIXEDX;
        mm[2] = All.DivBerrDensPdfDivBMin;
        mm[3] = All.DivBerrDensPdfDivBMax;

    }

    pdf2d( divbdp_x, divbdp_y, NULL, N, "DivBerrDensPdf", flag, mm );

    myfree( divbdp_x );
    myfree( divbdp_y );
    put_end();
}
#endif
