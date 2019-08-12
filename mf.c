#include "allvars.h"

void mass_function() {

    int g, *num, i, *cum_num, N1, N2, N;
    double m_min, m_max, M, dlogm,
           gm_min, gm_max, m_split, t;
    char buf[100];
    FILE *fd;

    writelog( "compute mass function ...\n" );

    gm_min = DBL_MAX;
    gm_max = -DBL_MAX;

    for ( g=0; g<Ngroups; g++ ) {
        vmin2( gm_min, Gprops[g].mass );
        vmax2( gm_max, Gprops[g].mass );
    }

    writelog( "gm_min: %g, gm_max: %g\n", gm_min, gm_max );

    m_min = gm_min;
    if ( All.MFMmin != 0 )
        m_min = All.MFMmin;

    m_max = gm_max;
    if ( All.MFMmax != 0 )
        m_max = All.MFMmax;

    writelog( "m_min: %g, m_max: %g\n", m_min, m_max );

    if ( m_min == m_max ) {
        endruns( "some error appear!" );
    }

    m_split = All.MFMSplit;

    t = pow( 10, (int)(log10( m_split )) );
    if ( t != m_split ) {
        writelog( "reset `m_split` from %.3e to %.3e\n", m_split,
             t );
        m_split = t;
    }

    N1 = All.MFBins;

    if ( m_max > m_split ) {
       // dlogm = log10( m_split/m_min ) / (N1-1);
        dlogm = log10( m_split/m_min ) / N1;
        N2 = (int)(log10(m_max)) - (int)(log10(m_split)) + 1;
    }
    else {
        N2 = 0;
        //dlogm = log10( m_max/m_min ) / (N1-1);
        dlogm = log10( m_split/m_min ) / N1;
    }

    writelog( "N1: %i, N2: %i\n", N1, N2 );

    N = N1+N2;

    mymalloc2( num, sizeof(int) * (N) );
    mymalloc2( cum_num, sizeof(int) * (N) );

    for ( g=0; g<Ngroups; g++ ) {

        M = Gprops[g].mass;

        if ( M < m_min || M > m_max )
            continue;

        if ( M >= m_split )
            i = N1 + log10( M / m_split );
        else
            i = log10(M/m_min) / dlogm;

        num[i] ++;

    }

    cum_num[N-1] = num[N-1];

    for ( i=N-2; i>=0; i-- ) {
        cum_num[i] = cum_num[i+1] + num[i];
    }

    sprintf( buf, "%s/MF/", OutputDir );
    create_dir( buf );
    sprintf( buf, "%s/MF/MF_%03i.dat", OutputDir, SnapIndex );

    fd = fopen( buf, "w" );

    for ( i=0; i<N1; i++ ) {

        fprintf( fd, "%e %g %g %i\n", pow( 10, i*dlogm + log10(m_min) ) * 1e10,
                num[i] / CUBE( BoxSize / All.MpcFlag ) / dlogm,
                cum_num[i] / CUBE( BoxSize / All.MpcFlag ),
                num[i]
                );

    }

    for ( i=0; i<N2; i++ ) {

        fprintf( fd, "%e %g %g %i\n", pow( 10, i ) * m_split * 1e10,
                num[N1+i] / CUBE( BoxSize / All.MpcFlag ),
                cum_num[N1+i] / CUBE( BoxSize / All.MpcFlag ),
                num[N1+i]
                );

    }

    fclose( fd );

    myfree( num );
    myfree( cum_num );

    sprintf( buf, "%s/MF/PS_%03i.dat", OutputDir, SnapIndex );

    fd = fopen( buf, "w" );

    for( M=m_min; M<=m_max; M*=1.1 )
        fprintf( fd, "%g %g\n", M*1e10, PS_dndM( Time, M ) * CUBE( All.MpcFlag ) * M );


    fclose( fd );

    writelog( "compute mass function ... done.\n" );
    put_sep;

}
