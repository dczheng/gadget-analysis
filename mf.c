#include "allvars.h"

void mass_function() {

    int g, *num, i, i_m_min, i_m_max;
    double m_min, m_max;
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

    i_m_min = floor( log10(m_min) );
    i_m_max = floor( log10(m_max) );

    writelog( "i_m_min: %i, i_m_max: %i\n", i_m_min, i_m_max );

    if ( i_m_min == i_m_max ) {
        endrun( "some error appear!" );
    }

    mymalloc2( num, sizeof(int) * ( i_m_max - i_m_min + 1 ) );

    for ( g=0; g<Ngroups; g++ ) {
        num[ ( int ) ( floor( log10( Gprops[g].mass * 1e10 ) ) - i_m_min ) ] ++;
    }

    for ( i=i_m_max-1; i>=i_m_min; i-- ) {
        num[i-i_m_min] += num[i-i_m_min+1];
    }

    /*
    for ( i=m_min; i<=m_max; i++ ) {
        writelog( "%i\n", num[i] );
    }
    */

    create_dir( "./mf" );
    sprintf( buf, "./mf/%s_mf_%.2f.dat", All.FilePrefix, All.RedShift );

    fd = fopen( buf, "w" );

    for ( i=i_m_min; i<=i_m_max; i++ ) {

        fprintf( fd, "%i %g\n", i, num[i-i_m_min] / CUBE( All.BoxSize / All.MpcFlag ) );

    }

    fclose( fd );

    myfree( num );

    writelog( "compute mass function ... done.\n" );
    put_block_line;
}
