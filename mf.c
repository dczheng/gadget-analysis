#include "allvars.h"

void mass_function() {
#ifdef MF

    int g, *num, i, *cum_num, N;
    double m_min, m_max, M, dlogm;
    FILE *fd;

    writelog( "compute mass function ...\n" );

    m_min = DBL_MAX;
    m_max = -DBL_MAX;

    for ( g=0; g<Ngroup; g++ ) {
        vmin2( m_min, Gprops[g].mass );
        vmax2( m_max, Gprops[g].mass );
    }


    if ( All.MFMmin != 0 )
        m_min = All.MFMmin;

    if ( All.MFMmax != 0 )
        m_max = All.MFMmax;

    writelog( "m_min: %g, m_max: %g\n", m_min, m_max );

    if ( m_min == m_max ) {
        endruns( "some error appear!" );
    }

    N = All.MFBins;
    dlogm = log10( m_max/m_min ) / (N-1);

    mymalloc2( num, sizeof(int) * (N) );
    mymalloc2( cum_num, sizeof(int) * (N) );

    for ( g=0; g<Ngroup; g++ ) {

        M = Gprops[g].mass;
        i = floor(log10( M/m_min ) / dlogm);
        if ( i<0 || i>N-1 )
            continue;
        num[i] ++;

    }

    cum_num[N-1] = num[N-1];

    for ( i=N-2; i>=0; i-- ) {
        cum_num[i] = cum_num[i+1] + num[i];
    }

    create_dir( "%s/MF/", OutputDir );
    fd = myfopen( "w", "%s/MF/MF_%03i.dat", OutputDir, SnapIndex );

    for ( i=0; i<N; i++ ) {

        fprintf( fd, "%e %g %g %i\n", m_min*pow(10, i*dlogm),
                num[i] / CUBE( BoxSize / All.MpcFlag ) / dlogm,
                cum_num[i] / CUBE( BoxSize / All.MpcFlag ),
                num[i]
                );

    }

    fclose( fd );
    myfree( num );
    myfree( cum_num );

    fd = myfopen( "w", "%s/MF/PS_%03i.dat", OutputDir, SnapIndex );
    for( M=m_min; M<=m_max; M*=1.1 )
        fprintf( fd, "%g %g\n", M,  PS_dndM( Time, M ) * CUBE( All.MpcFlag ) * M );
    fclose( fd );

    writelog( "compute mass function ... done.\n" );

#endif
}
