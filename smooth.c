#include "allvars.h"

//#define SMOOTH_DEBUG
//#define SMOOTH_DEBUG2
#define HIGH_PRECISION

void smooth() {

#if defined(BSMOOTH) || defined(CRESMOOTH) || defined (MACHSMOOTH)
    int ngbnum, k, n;
    long i, j;
    double h_i, h2_i, h_j, h2_j, r, cm[3], r2, t, hinv_i, hinv3_i, hinv4_i, u,
           wk, dwk, m_j, hinv_j, hinv3_j, hinv4_j;
#ifdef HIGH_PRECISION
    long double B[3], densitynorm, fac;
#else
     double B[3], densitynorm, fac;
#endif
#ifdef SMOOTH_DEBUG2
     long debug_i=99;
#endif

#ifdef BSMOOTH
    double *SmoothedB;
    mymalloc2( SmoothedB, sizeof(double) * N_Gas * 3 );
#endif

#ifdef MACHSMOOTH
    double *SmoothedM;
    mymalloc2( SmoothedM, sizeof(double) * N_Gas );
#endif

#ifdef CRESMOOTH
    double *SmoothedCRE_e, *SmoothedCRE_n;
    mymalloc2( SmoothedCRE_e, sizeof(double) * N_Gas );
    mymalloc2( SmoothedCRE_n, sizeof(double) * N_Gas );
#endif

    mytimer_start();
    put_header( "smooth" );
    writelog(
#ifdef BSMOOTH
    "B smooth\n"
#endif
#ifdef CRESMOOTH
    "CRE smooth\n"
#endif
#ifdef MACHSMOOTH
    "MACH smooth\n"
#endif
    );

#ifdef HIGH_PRECISION
    writelog( "long double: %li, double: %li\n",
        sizeof(long double), sizeof(double) );
#endif

    (void)h2_i;

    mymalloc1( Ngblist, NumPart * sizeof(long) );

    for( i=0; i<N_Gas; i++ ) {

#ifdef SMOOTH_DEBUG2
        i = debug_i;
        if (ThisTask_Local)
            continue;
#else
        if ( i % NTask_Local != ThisTask_Local )
            continue;
#endif

        h_i = SphP[i].Hsml;
        h2_i = h_i*h_i;
        kernel_hinv( h_i, &hinv_i, &hinv3_i, &hinv4_i );
        densitynorm = 0;
        for( k=0; k<3; k++ ) {
            cm[k] = P[i].Pos[k];
            B[k] = 0;
        }

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
                     B[k] += fac * SphP[j].B[k];
#endif
                }
#ifdef CRESMOOTH
                SmoothedCRE_n[i]      += fac * SphP[j].CRE_n;
                SmoothedCRE_e[i]      += fac * SphP[j].CRE_e;
#endif
#ifdef MACHSMOOTH
                SmoothedM[i] += fac * SphP[j].MachNumber;
#endif

                densitynorm += fac; 
            }
#ifdef SMOOTH_DEBUG2
            printf( "(%g %g %g): %g\n",
            P[j].Pos[0],
            P[j].Pos[1],
            P[j].Pos[2],
            get_B(j)
            );
#endif

        } // for n

        if ( densitynorm > 0 ) {
            for( k=0; k<3; k++ ) {
#ifdef BSMOOTH
                SmoothedB[i*3+k] = B[k] / densitynorm;
#endif
            }
#ifdef CRESMOOTH
            SmoothedCRE_e[i]      /= densitynorm; 
            SmoothedCRE_n[i]      /= densitynorm; 
#endif
#ifdef MACHSMOOTH
            SmoothedM[i] /= densitynorm;
#endif
        }
#ifdef SMOOTH_DEBUG2
        break;
#endif
    } // for i

    do_sync_local( "smooth" );

    for ( i=0; i<N_Gas; i++ ) {
        if ( i % NTask_Local != ThisTask_Local )
            continue;
        for( k=0; k<3; k++ ) {
#ifdef BSMOOTH
            SphP[i].B[k] = SmoothedB[i*3+k];
#endif
        }

#ifdef CRESMOOTH
        SphP[i].CRE_n      =  SmoothedCRE_n[i];
        SphP[i].CRE_e      =  SmoothedCRE_e[i];
#endif
#ifdef MACHSMOOTH
        SphP[i].MachNumber =  SmoothedM[i];
#endif
    }

    myfree( Ngblist );
    mytimer_end();

#ifdef BSMOOTH
    myfree( SmoothedB );
#endif
#ifdef MACHSMOOTH
    myfree( SmoothedM )
#endif
#ifdef CRESMOOTH
    myfree( SmoothedCRE_n );
    myfree( SmoothedCRE_e );
#endif

    do_sync_local( "smooth" );

#ifdef SMOOTH_DEBUG
    if ( !ThisTask_Local ) {
#ifdef SMOOTH_DEBUG2
        printf( "[%li] %g\n", debug_i, get_B(debug_i) );
#else
        for( i=0; i<100; i++ )
            printf( "[%li] %g\n", i, get_B(i) );
#endif
        endruns( "smooth-test" );
    }
#endif
    put_end();

#endif

}

void smooth2() {

#if defined(RADSMOOTH)
    int ngbnum, k, n;
    long i, j;
    double h_i, h2_i, h_j, h2_j, r, cm[3], r2, t, hinv_i, hinv3_i, hinv4_i, u,
           wk, dwk, m_j, hinv_j, hinv3_j, hinv4_j, densitynorm, fac;

#ifdef RADSMOOTH
    double *Smoothedrad;
    mymalloc2( Smoothedrad, sizeof(double) * N_Gas * All.FreqN );
#endif

    mytimer_start();
    put_header( "smooth2" );
    writelog(
#ifdef RADSMOOTH
    "rad smooth\n"
#endif
    );

    (void)h2_i;

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
        }

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
#ifdef RADSMOOTH
                for( k=0; k<All.FreqN; k++ )
                    Smoothedrad[i*All.FreqN+k] += fac * PartRad[i*All.FreqN+k];
#endif

                densitynorm += fac; 
            }

        } // for n

        if ( densitynorm > 0 ) {
            for( k=0; k<All.FreqN; k++ ) {
                Smoothedrad[i*All.FreqN+k] /= densitynorm;
            }
        }
    } // for i

    do_sync_local( "smooth2" );

    for ( i=0; i<N_Gas; i++ ) {
        if ( i % NTask_Local != ThisTask_Local )
            continue;
        for( k=0; k<All.FreqN; k++ ) {
            PartRad[i*All.FreqN+k] = Smoothedrad[i*All.FreqN+k];
        }

    }

    myfree( Ngblist );
    mytimer_end();

#ifdef RADSMOOTH
    myfree( Smoothedrad );
#endif

    do_sync_local( "smooth2" );

    put_end();

#endif

}
