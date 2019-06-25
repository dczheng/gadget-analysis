#include "allvars.h"

void part_info() {

    long i;
    double rho_max, rho_min;

    writelog( "Particle info\n" );

    rho_max = -DBL_MAX;
    rho_min = DBL_MAX;

    for( i=0; i<N_Gas; i++ ) {
        vmax2( rho_max, SphP[i].Density );
        vmin2( rho_min, SphP[i].Density, 1 );
    }

    writelog( "[Baryon density] min: %g [%g gcm-3], max: %g [%g gcm-3]\n"
              "[Proton number density], min: %g cm-3, max: %g cm-3\n"
              "RhoBar: %g gcm-3, RhoProBar: %g cm-3\n",
            rho_min,
            rho_min * g2c.g / CUBE( g2c.cm ),
            rho_max,
            rho_max * g2c.g / CUBE( g2c.cm ),
            rho_min * g2c.g / CUBE( g2c.cm ) / cuc.m_p,
            rho_max * g2c.g / CUBE( g2c.cm ) / cuc.m_p,
            All.RhoBaryon * g2c.g / CUBE( g2c.cm ),
            All.RhoBaryon * g2c.g / CUBE( g2c.cm ) / cuc.m_p
            );

    put_sep;
}

