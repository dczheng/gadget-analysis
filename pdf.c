#include"allvars.h"

void dens_pdf() {
#ifdef DENSPDF

    double *num, *num_warm_hot, dlogDens, Densmin, Densmax, *num_diffuse_condense, *num_hot;
    int N, i;
    long p;
    FILE *fd;
    
    N = All.DensPdfN;

    get_gas_density_min_max( Densmin, Densmax );

    Densmin /= RhoBaryon;
    Densmax /= RhoBaryon;

    if ( All.DensPdfMin > 0 )
        Densmin = All.DensPdfMin;
    if ( All.DensPdfMax > 0 )
        Densmax = All.DensPdfMax;

    writelog( "Den, min: %g, max: %g\n", Densmin, Densmax );

    dlogDens = log10( Densmax / Densmin ) / N;
    mymalloc2( num, sizeof(double) * N );
    mymalloc2( num_warm_hot, sizeof(double) * N );
    mymalloc2( num_hot, sizeof(double) * N );
    mymalloc2( num_diffuse_condense, sizeof(double) * N );

    for( p=0; p<N_Gas; p++ ) {
        i = log10( SphP[p].Density / Time3 / (Densmin * RhoBaryon) ) / dlogDens;
        if ( i >= N || i<0 )
            continue;

        num[i] ++;

        if ( SphP[p].Temp <= 1e5 )
            num_diffuse_condense[i] ++;

        if ( SphP[p].Temp > 1e5 && SphP[p].Temp < 1e7 )
            num_warm_hot[i] ++;

        if ( SphP[p].Temp >= 1e7 )
            num_hot[i] ++;
    }

    for( i=0; i<N; i++ ) {
        num[i] /= dlogDens;
        num_diffuse_condense[i] /= dlogDens;
        num_hot[i] /= dlogDens;
    }

    create_dir( "%s/DensPdf", OutputDir );
    fd = myfopen( "w", "%s/DensPdf/DensPdf_%03i.dat", OutputDir, SnapIndex );

    fprintf( fd, "%g 0 0 0 0\n",  Redshift );
    for( i=0; i<N; i++ )
        fprintf( fd, "%g %g %g %g %g\n", Densmin*pow( 10, i*dlogDens ),
        num[i], 
        num_diffuse_condense[i],
        num_warm_hot[i],
        num_hot[i]
        );

    fclose( fd );
    myfree( num );
    myfree( num_warm_hot );
    myfree( num_diffuse_condense );
    myfree( num_hot );

#endif
}

void cree_pdf() {
#ifdef CREEPDF

    double *num, dloge, emin, emax, e;
    int N, i;
    long p;
    FILE *fd;
    
    N = All.CREePdfN;

    emin = DBL_MAX;
    emax = -emin;

    for( p=0; p<N_Gas; p++ ) {
        e = SphP[p].CRE_e / SphP[p].u;
        vmin20( emin, e );
        vmax2( emax, e );
    }

    if ( All.CREePdfMin > 0 )
        emin = All.CREePdfMin;
    if ( All.CREePdfMax > 0 )
        emax = All.CREePdfMax;

    writelog( "CREe, min: %g, max: %g\n", emin, emax );

    dloge = log10( emax / emin ) / (N-1);
    mymalloc2( num, sizeof(double) * N );

    for( p=0; p<N_Gas; p++ ) {
        e = SphP[p].CRE_e / SphP[p].u;
        i = floor(log10( e/emin ) / dloge);
        if ( i > N-1 || i < 0 )
            continue;
        num[i] ++;
    }

    for( i=0; i<N; i++ )
        num[i] = num[i] / ((double)N_Gas) / dloge;

    create_dir( "%s/CREePdf", OutputDir );
    fd = myfopen( "w", "%s/CREePdf/CREePdf_%03i.dat", OutputDir, SnapIndex );

    fprintf( fd, "%g %g\n", Redshift, Redshift );
    for( i=0; i<N; i++ )
        fprintf( fd, "%g %g\n", emin*pow( 10, i*dloge ), num[i] );

    fclose( fd );
    myfree( num );

#endif
}

void T_pdf() {
#ifdef TPDF

    double *num, dlogT, Tmin, Tmax, *num_fof, s, s_fof, *num4, *num3, s4, s3;
    int N, i, j;
    long p;
    FILE *fd;
    struct group_properties *g;
    
    N = All.TPdfN;

    put_header( "T Pdf" );
    get_gas_temp_min_max( Tmin, Tmax );

    //printf( "%g %g\n", All.TPdfMin, All.TPdfMax );

    if ( All.TPdfMin > 0 )
        Tmin = All.TPdfMin;

    if ( All.TPdfMax > 0 )
        Tmax = All.TPdfMax;

    writelog( "T, min: %g, max: %g\n", Tmin, Tmax );

    dlogT = log10( Tmax / Tmin ) / N;

    mymalloc2( num, sizeof(double) * N );
    mymalloc2( num_fof, sizeof(double) * N );
    mymalloc2( num3, sizeof(double) * N );
    mymalloc2( num4, sizeof(double) * N );

    for(p=0; p<NumPart; p++)
        P[p].flag = 0;

    for(i=0; i<Ngroup; i++) {
        g = &Gprops[i];
        p = g->Head;
        while( p>=0 ) {

            j = floor(log10( SphP[p].Temp / Tmin ) / dlogT);

            if ( P[p].Type ||
                j > N-1 ||
                j < 0 ) {
                p = FoFNext[p];
                continue;
            }

            if ( g->mass > 1e4 ) {
                num4[j] ++;
                s4++;
            }

            if ( g->mass <= 1e4 && g->mass > 1e3 ) {
                num3[j] ++;
                s3++;
            }

            p = FoFNext[p];
        }
    }

    for(i=0; i<Ngroup; i++) {
        g = &Gprops[i];
        p = g->Head;
        if ( g->mass < 10 )
            continue;
        while( p>=0 ) {
            P[p].flag = 1;
            p = FoFNext[p];
        }
    }


    s = s_fof = 0;
    for( p=0; p<N_Gas; p++ ) {
        j = floor(log10( SphP[p].Temp / Tmin ) / dlogT);
        if ( j > N-1 || j < 0 )
            continue;
        if ( P[p].flag ) {
            num_fof[j] ++;
            s_fof ++;
        }
        else {
            num[j] ++;
            s++;
        }

    }

    create_dir( "%s/TPdf", OutputDir );
    fd = myfopen( "w", "%s/TPdf/TPdf_%03i_num.dat", OutputDir, SnapIndex );

    fprintf( fd, "%g 0 0  0 0 \n", Redshift );
    for( i=0; i<N; i++ )
        fprintf( fd, "%g %g %g %g %g\n", Tmin*pow( 10, i*dlogT ), num_fof[i], num[i],
            num4[i], num3[i] );

    fclose( fd );

    for( i=0; i<N; i++ ) {
        num[i] /= dlogT*s;
        num_fof[i] /= dlogT*s_fof;
        num4[i] /= dlogT*s4;
        num3[i] /= dlogT*s3;
    }

    fd = myfopen( "w", "%s/TPdf/TPdf_%03i.dat", OutputDir, SnapIndex );

    fprintf( fd, "%g 0 0  0 0 \n", Redshift );
    for( i=0; i<N; i++ )
        fprintf( fd, "%g %g %g %g %g\n", Tmin*pow( 10, i*dlogT ), num_fof[i], num[i],
            num4[i], num3[i] );

    fclose( fd );
    myfree( num );
    myfree( num_fof );
    myfree( num4 );
    myfree( num3 );
    put_end();
#endif

}

void B_Pdf() {
#ifdef BPDF
    long i;
    int j, N;
    double Bmin, Bmax, *n, dlogB, B, s, *nn;
    FILE *fd;

//#define B_PDF_DEBUG

#ifdef B_PDF_DEBUG
    int nnn=0;
#endif
    put_header( "B pdf" );
    get_B_min_max( Bmin, Bmax );
    N = All.BPdfBins;

    if ( All.BPdfBMin > 0 )
        Bmin = All.BPdfBMin;

    if ( All.BPdfBMax > 0 )
        Bmax = All.BPdfBMax;

    dlogB = log10(Bmax/Bmin) / N;

    mymalloc2( n, sizeof(double) * N );
    mymalloc2( nn, sizeof(double) * N );

    for( i=0; i<N_Gas; i++ ) {
        B = get_B( i );
        if ( B == 0 )
            continue;

        j = log10( B/Bmin ) / dlogB;
        if ( j<0 || j>N-1 )
            continue;
        n[j] ++;
#ifdef B_PDF_DEBUG
        nnn ++; 
#endif
    }
#ifdef B_PDF_DEBUG
    printf( "task: %i, nnn: %i\n", ThisTask, nnn );
#endif

    for( i=0, s=0; i<N; i++ ) {
        s += n[i];
    }

    for( i=0; i<N; i++ )
        nn[i] = n[i] / s / dlogB;

    create_dir( "%s/BPdf", OutputDir );
    fd = myfopen( "w", "%s/BPdf/BPdf_%03i.dat", OutputDir, SnapIndex );

    fprintf( fd, "%g 0 0\n",  Redshift );
    for( i=0; i<N; i++ )
        fprintf( fd, "%g %g %g\n", Bmin * pow(10, i*dlogB) * 1e6, nn[i], n[i] );

    fclose( fd );
    myfree( n );
    myfree( nn );

#ifdef B_PDF_PRINT_LARGE_B 
    int task;
    for( task=0; task<NTask; task++ ) {
        if ( task == ThisTask ) {
            for( j=10,N=0; j<1000; j*=10 ) {
                for( i=0; i<N_Gas; i++ ) {
                    B = get_B( i ) * 1e6;
                    if ( B>j && B<j*10 ) {
                        printf( "[%03i, %05i], B: %e, density: %e, temp: %e, sfr: %g\n",
                        task, N, B, 
                        SphP[i].Density / Time3 / RhoBaryon,
                        SphP[i].Temp,
                        SphP[i].sfr
                        );
                        N++;
                    }
                }
            }
        }
        MPI_Barrier( MpiComm_Master );
    }
#endif

/*
    for( i=0, N=0; i<N_Gas; i++ ) {
        if ( SphP[i].sfr > 0 ) {
            B = get_B( i ) * 1e6;
            writelog( "[%05i], B: %e, density: %e, temp: %e, sfr: %g\n",
            N, B, 
            SphP[i].Density / Time3 / RhoBaryon,
            SphP[i].Temp,
            SphP[i].sfr
            );
            N++;
        }
    }
*/

#endif
}

void DivB_Err_Pdf() {
#ifdef DIVBERRPDF
    long i, idxmin, idxmax;
    int j, N;
    double errmin, errmax, *n, dlogerr, err, B, errmean, vsum, v;
    FILE *fd;

    put_header( "DivB Error pdf" );
    errmin = 1e100;
    errmax = -errmin;
    for( i=0; i<N_Gas; i++ )  {
        B = get_B(i);
        if ( B<=0 )
            continue;
        err = SphP[i].divB * SphP[i].Hsml / B;
        err = fabs(err);
        if ( err>errmax ) {
            errmax = err;
            idxmax = i;
        }
        if ( err<errmin ) {
            errmin = err;
            idxmin = i;
        }
    }

    N = All.DivBErrPdfBins;

    if ( All.DivBErrPdfMin > 0 )
        errmin = All.DivBErrPdfMin;

    if ( All.DivBErrPdfMax > 0 )
        errmax = All.DivBErrPdfMax;
    //writelog( "N: %i, DivBErrMin: %g, divBErrMax: %g\n", N, errmin, errmax );

    dlogerr = log10(errmax/errmin) / N;

    mymalloc2( n, sizeof(double) * N );

    errmean = vsum = 0;
    for( i=0; i<N_Gas; i++ ) {
        B = get_B( i );
        if ( B <= 0 )
            continue;

        err = SphP[i].divB * SphP[i].Hsml / B;
        err = fabs(err);
        v = P[i].Mass / SphP[i].Density;
        errmean += err * v;
        vsum += v;
        j = log10( err/errmin ) / dlogerr;
        if ( j<0 || j>N-1 )
            continue;
        n[j] ++;
    }

    errmean /= vsum;
    writelog( "errmin: %g [%g %g %g %g]\nerrmax: %g [%g %g %g %g]\nerrmean: %g\n",
        errmin, SphP[idxmin].divB, SphP[idxmin].Hsml, get_B( idxmin ),
        SphP[idxmin].Density / Time3 / RhoBaryon,
        errmax, SphP[idxmax].divB, SphP[idxmin].Hsml, get_B( idxmax ),
        SphP[idxmax].Density / Time3 / RhoBaryon,
        errmean );

    for( i=0; i<N; i++ )
        n[i] = n[i] / dlogerr;

    create_dir( "%s/DivBErrPdf", OutputDir );
    fd = myfopen( "w", "%s/DivBErrPdf/DivBErrPdf_%03i.dat", OutputDir, SnapIndex );

    fprintf( fd, "%g \n",  errmean );
    for( i=0; i<N; i++ )
        fprintf( fd, "%g %g\n", errmin * pow(10, i*dlogerr), n[i] );

    fclose( fd );
    myfree( n );

#endif
}
void crep_pdf() {

#ifdef CREPPDF
    double *num, dlogp, pmin, pmax, pp;
    int N, i;
    long p;
    FILE *fd;
    
    N = All.CREPPdfN;

    pmin = DBL_MAX;
    pmax = -pmin;

    for( p=0; p<N_Gas; p++ ) {
        pp = SphP[p].CRE_P / (SphP[p].CRE_P + SphP[p].CR_P0 +
                get_pressure(p) );
        vmin20( pmin, pp );
        vmax2( pmax, pp );
    }

    if ( All.CREPPdfMin > 0 )
        pmin = All.CREPPdfMin;
    if ( All.CREPPdfMax > 0 )
        pmax = All.CREPPdfMax;
    writelog( "CREP, min: %g, max: %g\n", pmin, pmax );

    dlogp = log10( pmax / pmin ) / (N-1);
    mymalloc2( num, sizeof(double) * N );

    for( p=0; p<N_Gas; p++ ) {
        pp = SphP[p].CRE_P / (SphP[p].CRE_P + SphP[p].CR_P0 +
                get_pressure(p) );
        i = floor(log10( pp/pmin ) / dlogp);
        if ( i > N-1 || i < 0 )
            continue;
        num[i] ++;
    }

    for( i=0; i<N; i++ )
        num[i] = num[i] / ((double)N_Gas) / dlogp;

    create_dir( "%s/CREPPdf", OutputDir );
    fd = myfopen( "w", "%s/CREPPdf/CREPPdf_%03i.dat", OutputDir, SnapIndex );

    fprintf( fd, "%g %g\n", Redshift, Redshift );
    for( i=0; i<N; i++ )
        fprintf( fd, "%g %g\n", pmin*pow( 10, i*dlogp ), num[i] );

    fclose( fd );
    myfree( num );
#endif

}

void crp_pdf() {

#ifdef CRPPDF
    double *num, dlogp, pmin, pmax, pp;
    int N, i;
    long p;
    FILE *fd;
    
    N = All.CREPPdfN;

    pmin = DBL_MAX;
    pmax = -pmin;

    for( p=0; p<N_Gas; p++ ) {
        pp = SphP[p].CR_P0 / (SphP[p].CRE_P + SphP[p].CR_P0 +
                get_pressure(p) );
        vmin20( pmin, pp );
        vmax2( pmax, pp );
    }

    if ( All.CRPPdfMin > 0 )
        pmin = All.CRPPdfMin;
    if ( All.CRPPdfMax > 0 )
        pmax = All.CRPPdfMax;
    writelog( "CREP, min: %g, max: %g\n", pmin, pmax );

    dlogp = log10( pmax / pmin ) / (N-1);
    mymalloc2( num, sizeof(double) * N );

    for( p=0; p<N_Gas; p++ ) {
        pp = SphP[p].CR_P0 / (SphP[p].CRE_P + SphP[p].CR_P0 +
                get_pressure(p) );
        i = floor(log10( pp/pmin ) / dlogp);
        if ( i > N-1 || i < 0 )
            continue;
        num[i] ++;
    }

    for( i=0; i<N; i++ )
        num[i] = num[i] / ((double)N_Gas) / dlogp;

    create_dir( "%s/CRPPdf", OutputDir );
    fd = myfopen( "w", "%s/CRPPdf/CRPPdf_%03i.dat", OutputDir, SnapIndex );

    fprintf( fd, "%g %g\n", Redshift, Redshift );
    for( i=0; i<N; i++ )
        fprintf( fd, "%g %g\n", pmin*pow( 10, i*dlogp ), num[i] );

    fclose( fd );
    myfree( num );
#endif

}

void mach_pdf() {

#ifdef MACHPDF
    double *num, dlogm, mmin, mmax, m, v, vsum;
    int N, i;
    long p;
    FILE *fd;
    
    N = All.MachPdfN;

    mmin = DBL_MAX;
    mmax = -mmin;

    for( p=0; p<N_Gas; p++ ) {
        m = SphP[p].MachNumber;
        vmin20( mmin, m );
        vmax2( mmax, m );
    }

    if ( All.MachPdfMin > 0 )
        mmin = All.MachPdfMin;
    if ( All.MachPdfMax > 0 )
        mmax = All.MachPdfMax;
    writelog( "Mach, min: %g, max: %g\n", mmin, mmax );

#define USE_LOG
#ifdef USE_LOG
    dlogm = log( mmax / mmin ) / (N-1);
#else
    dlogm = ( mmax - mmin ) / (N-1);
#endif
    mymalloc2( num, sizeof(double) * N );

    for( p=0,vsum; p<N_Gas; p++ ) {
        m = SphP[p].MachNumber;
        v = get_V(p);
        //v = 4.0 / 3.0 * PI * P[p].Mass / SphP[p].Density;
        //v = CUBE( SphP[p].Hsml );
#ifdef USE_LOG
        i = floor(log( m/mmin ) / dlogm);
#else
        i = floor( (m-mmin) / dlogm);
#endif
        if ( i > N-1 || i < 0 )
            continue;
        num[i] += v;
        vsum += v;
    }

    for( i=0; i<N; i++ ) {
        //num[i] = num[i] / ((double)N_Gas);
        //num[i] = num[i] / dlogm;
        num[i] = num[i] / vsum / dlogm;
        //num[i] = num[i] / CUBE(BoxSize);
    }

    create_dir( "%s/MachPdf", OutputDir );
    fd = myfopen( "w", "%s/MachPdf/MachPdf_%03i.dat", OutputDir, SnapIndex );

    fprintf( fd, "%g %g\n", Redshift, Redshift );
    for( i=0; i<N; i++ )
#ifdef USE_LOG
        fprintf( fd, "%g %g\n", mmin*exp( i*dlogm ), num[i] );
#else
        fprintf( fd, "%g %g\n", mmin +  i*dlogm, num[i] );
#endif

    fclose( fd );
    myfree( num );
#endif

}
