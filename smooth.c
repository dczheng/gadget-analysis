#include "allvars.h"


void smooth() {
#ifdef SMOOTH

//#define SMOOTH_DEBUG
    int ngbnum, k, n;
    long i, j;
    double h_i, h2_i, h_j, h2_j, r, cm[3], r2, t, hinv_i, hinv3_i, hinv4_i, u,
           wk, dwk, m_j, densitynorm, hinv_j, hinv3_j, hinv4_j;

    (void)h2_i;

    mytimer_start();
    put_header( "smooth" );
#ifdef BSMOOTH
    writelog( "B smooth\n" );
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

                for( k=0; k<3; k++ ) {
#ifdef BSMOOTH
                     SphP[i].SmoothB[k] += m_j * wk * SphP[j].B[k];
#endif
                }

                densitynorm += m_j * wk; 
            }

        }
#ifdef BSMOOTH
        if ( densitynorm > 0 ) {
            for( k=0; k<3; k++ )
                SphP[i].SmoothB[k] /= densitynorm;
                
        }
#endif
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
    }

    myfree( Ngblist );

    mytimer_end();

#ifdef SMOOTH_DEBUG
    endruns( "smooth-test" );
#endif
    put_end();

#endif
}
