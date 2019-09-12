#include "allvars.h"

int smooth_flag() {

    if ( All.BSmooth )
        return 1;
    return 0;
}

void smooth() {

//#define SMOOTH_DEBUG
    int ngbnum, k, n;
    long i, j;
    double h_i, h2_i, h_j, h2_j, r, cm[3], r2, t, hinv_i, hinv3_i, hinv4_i, u,
           wk, dwk, m_j, densitynorm, hinv_j, hinv3_j, hinv4_j;

    (void)h2_i;
    

    if ( !smooth_flag() )
        return;

    mytimer_start();
    writelog( "smooth ...\n" );
    mymalloc1( Ngblist, NumPart * sizeof(long) );

    for( i=0; i<N_Gas; i++ ) {
        h_i = SphP[i].Hsml;
        h2_i = h_i*h_i;
        kernel_hinv( h_i, &hinv_i, &hinv3_i, &hinv4_i );
        densitynorm = 0;
        for( k=0; k<3; k++ ) {
            cm[k] = P[i].Pos[k];
            if ( All.BSmooth ) {
                SphP[i].SmoothB[k] = 0;
            }
        }

        ngbnum = ngb( cm, h_i );
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
                    if ( All.BSmooth )
                        SphP[i].SmoothB[k] += m_j * wk * SphP[j].B[k];
                }

                densitynorm += m_j * wk; 
            }

        }
        if ( All.BSmooth && densitynorm > 0 ) {
            for( k=0; k<3; k++ )
                SphP[i].SmoothB[k] /= densitynorm;
                
        }
    }

    for ( i=0; i<N_Gas; i++ ) {
        if ( All.BSmooth )
            for( k=0; k<3; k++ )
                SphP[i].B[k] = SphP[i].SmoothB[k];
    }

    myfree( Ngblist );

    mytimer_end();

#ifdef SMOOTH_DEBUG
    endruns( "smooth-test" );
#endif
    put_sep0;
}
