#include "allvars.h"

void mass_function() {

    int g, *num, i, *cum_num, N;
    double m_min, m_max, M, dlogm;
    char buf[100];
    FILE *fd;
    writelog( "compute mass function ...\n" );

    m_min = DBL_MAX;
    m_max = -DBL_MAX;

    for ( g=0; g<Ngroups; g++ ) {
        m_min = vmin( m_min, Gprops[g].mass, 0 );
        m_max = vmax( m_max, Gprops[g].mass );
    }

    m_min *= 1e10;
    m_max *= 1e10;

    writelog( "m_min: %g, m_max: %g\n", m_min, m_max );

    if ( m_min == m_max ) {
        endrun( "some error appear!" );
    }

    N = All.MFBins;
    dlogm = log10( m_max/m_min ) / (N-1);

    mymalloc2( num, sizeof(int) * N );
    mymalloc2( cum_num, sizeof(int) * N );


    for ( g=0; g<Ngroups; g++ ) {
        M = Gprops[g].mass * 1e10;
        num[ (int)(log10(M/m_min) / dlogm) ] ++;
    }

    cum_num[N-1] = num[N-1];

    for ( i=-2; i>=0; i-- ) {
        cum_num[i] = cum_num[i+1] + num[i];
    }

    create_dir( "./MF" );
    sprintf( buf, "./MF/MF_%.2f.dat", All.RedShift );

    fd = fopen( buf, "w" );

    for ( i=0; i<N; i++ ) {

        fprintf( fd, "%e %g %g\n", pow( 10, i*dlogm + log10(m_min) ),
                num[i] / CUBE( All.BoxSize / All.MpcFlag ) / dlogm,
                cum_num[i] / CUBE( All.BoxSize / All.MpcFlag )
                );

    }

    fclose( fd );

    myfree( num );
    myfree( cum_num );

    sprintf( buf, "./MF/PS_%.2f.dat", All.RedShift );

    fd = fopen( buf, "w" );

    for( M=1; M<1e6; M*=1.1 )
        fprintf( fd, "%g %g\n", M * 1e10, PS_dndM( All.Time, M ) * CUBE( All.MpcFlag ) / 1e10 * (M*1e10) );


    fclose( fd );


    writelog( "compute mass function ... done.\n" );
    put_block_line;
}
