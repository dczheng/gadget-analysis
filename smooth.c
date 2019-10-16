#include "allvars.h"


void smooth() {

#if defined(BSMOOTH) || defined(CRESMOOTH)
    int ngbnum, k, n;
    long i, j;
    double h_i, h2_i, h_j, h2_j, r, cm[3], r2, t, hinv_i, hinv3_i, hinv4_i, u,
           wk, dwk, m_j, densitynorm, hinv_j, hinv3_j, hinv4_j, fac;

    (void)h2_i;

    mytimer_start();
    put_header( "smooth" );
#ifdef BSMOOTH
    writelog( "B smooth\n" );
#endif
#ifdef CRESMOOTH
    writelog( "CRE smooth\n" );
#endif
    mymalloc1( Ngblist, NumPart * sizeof(long) );

    for( i=0; i<N_Gas; i++ ) {

        if ( i % NTask_Local != ThisTask_Local )
            continue;

        h_i = SphP[i].Hsml;
        h2_i = h_i*h_i;
        kernel_hinv( h_i, &hinv_i, &hinv3_i, &hinv4_i );
        densitynorm = 0;
        for( k=0; k<3; k++ ) {
            cm[k] = P[i].Pos[k];
#ifdef BSMOOTH
            SphP[i].SmoothB[k] = 0;
#endif
        }
#ifdef CRESMOOTH
        SphP[i].SmoothCRE_C =
        SphP[i].SmoothCRE_Alpha =
        SphP[i].SmoothCRE_qmin =
        SphP[i].SmoothCRE_qmax =
        SphP[i].SmoothCRE_n =
        SphP[i].SmoothCRE_e = 0;
#endif

        ngbnum = ngb( cm, h_i, 0 );
        if ( ngbnum == 0 )
            continue;

        for( n=0; n<ngbnum; n++ ) {

            j = Ngblist[n];
            if ( P[j].Type != 0 )
                continue;

            for( k=0, r2=0; k<3; k++ ) {
                t = P[j].Pos[k] - cm[k];
                t = PERIODIC_HALF( t );
                r2 += t*t;
            }

            r = sqrt( r2 );

            h_j = SphP[j].Hsml; 
            h2_j = h_j * h_j;

            if ( r2<h2_j ) {

                kernel_hinv( h_j, &hinv_j, &hinv3_j, &hinv4_j );
                m_j = P[j].Mass;

                u = r * hinv_j;
                kernel_main( u, hinv3_j, hinv4_j, &wk, &dwk, -1 );
                wk /= SphP[j].Density;

                fac = m_j * wk;
                for( k=0; k<3; k++ ) {
#ifdef BSMOOTH
                     SphP[i].SmoothB[k] += fac * SphP[j].B[k];
#endif
                }

#ifdef CRESMOOTH
                SphP[i].SmoothCRE_C      += fac * SphP[j].CRE_C;
                SphP[i].SmoothCRE_Alpha  += fac * SphP[j].CRE_Alpha;
                SphP[i].SmoothCRE_qmin   += fac * SphP[j].CRE_qmin;
                SphP[i].SmoothCRE_qmax   += fac * SphP[j].CRE_qmax;
                SphP[i].SmoothCRE_n      += fac * SphP[j].CRE_n;
                SphP[i].SmoothCRE_e      += fac * SphP[j].CRE_e;
#endif
                densitynorm += fac; 
            }

        }
        if ( densitynorm > 0 ) {
#ifdef BSMOOTH
            for( k=0; k<3; k++ )
                SphP[i].SmoothB[k] /= densitynorm;
#endif
#ifdef CRESMOOTH
            SphP[i].CRE_C      /= densitynorm; 
            SphP[i].CRE_Alpha  /= densitynorm; 
            SphP[i].CRE_qmin   /= densitynorm; 
            SphP[i].CRE_qmax   /= densitynorm; 
            SphP[i].CRE_n      /= densitynorm; 
            SphP[i].CRE_e      /= densitynorm; 
#endif
        }
    }

    for ( i=0; i<N_Gas; i++ ) {
        if ( i % NTask_Local != ThisTask_Local )
            continue;
#ifdef BSMOOTH
            /*
        if ( i < 100 ) {
            printf( "B: " );
            for( k=0; k<3; k++ )
                printf( "%g ", SphP[i].B[k] * 1e6 );
            printf( "SmoothB: " );
            for( k=0; k<3; k++ )
                printf( "%g ", SphP[i].SmoothB[k] * 1e6 );
            printf( "\n" );
        }
            */
        for( k=0; k<3; k++ )
            SphP[i].B[k] = SphP[i].SmoothB[k];
#endif

#ifdef CRESMOOTH
        SphP[i].CRE_C      =  SphP[j].SmoothCRE_C;
        SphP[i].CRE_Alpha  =  SphP[j].SmoothCRE_Alpha;
        SphP[i].CRE_qmin   =  SphP[j].SmoothCRE_qmin;
        SphP[i].CRE_qmax   =  SphP[j].SmoothCRE_qmax;
        SphP[i].CRE_n      =  SphP[j].SmoothCRE_n;
        SphP[i].CRE_e      =  SphP[j].SmoothCRE_e;
#endif
    }

    myfree( Ngblist );

    mytimer_end();

#ifdef SMOOTH_DEBUG
    endruns( "smooth-test" );
#endif
    put_end();

#endif

}
