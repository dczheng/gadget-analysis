#include "allvars.h"

void group_gas_ratio() {

    long p;
    struct group_properties g;
    int g_index, i;
    FILE *fd;

    double r_gas, r_condensed, r_diffuse, r_hot, r_warmhot, m_tot, m_gas;

    create_dir( "%sGroupGasRatio", OutputDir );
    fd = myfopen( "w", "%s/GroupGasRatio/GroupGasRatio_%03i.dat", OutputDir, SnapIndex );

    for ( g_index=0; g_index<Ngroups; g_index++ ) {
        if ( !group_present( g_index ) )
            break;
        g = Gprops[g_index];

        for( i=0,m_tot=0; i<6; i++ )
            m_tot += g.mass_table[i];
        m_gas = g.mass_table[0];
        r_gas = m_gas / m_tot;

        r_condensed = r_diffuse = r_hot = r_warmhot = 0; 
        p = g.Head;
        while( p >= 0 ) {
            if ( P[p].Type != 0 ) {
                p = FoFNext[p];
                continue;
            }

            if ( SphP[p].Temp >= 1e7 )
                r_hot += P[p].Mass;

            if ( SphP[p].Temp < 1e7 && SphP[p].Temp >= 1e5 )
                r_warmhot += P[p].Mass;

            if ( SphP[p].Temp < 1e5 && SphP[p].Density / Time3 / RhoBaryon < 1e3 )
                r_diffuse += P[p].Mass;

            if ( SphP[p].Temp < 1e5 && SphP[p].Density / Time3 / RhoBaryon >= 1e3 )
                r_condensed += P[p].Mass;

            p = FoFNext[p];
        }

        r_hot /= m_gas;
        r_warmhot /= m_gas; 
        r_condensed /= m_gas;
        r_diffuse /= m_gas;

        fprintf( fd, "%i %g %g %g %g %g\n", g_index,
            r_gas, r_hot, r_warmhot, r_condensed, r_diffuse );

    }

    fclose( fd );
}


