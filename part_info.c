#include "allvars.h"

void part_info() {

    long i;
    double rho_max, rho_min;
#ifdef READMACH
    double Mmax, Mmin;
#endif

    writelog( "Particle info\n" );

    rho_max = -DBL_MAX;
    rho_min = DBL_MAX;

#ifdef READMACH
    Mmax = rho_max;
    Mmin = rho_min;
#endif
    for( i=0; i<N_Gas; i++ ) {
        vmax2( rho_max, SphP[i].Density );
        vmin20( rho_min, SphP[i].Density );
#ifdef READMACH
        vmax2( Mmax, SphP[i].MachNumber );
        vmin20( Mmin, SphP[i].MachNumber );
#endif
    }

    writelog( "[Baryon density]\nmin: %g [%g gcm-3]\nmax: %g [%g gcm-3]\n"
              "[Proton number density]\nmin: %g cm-3\nmax: %g cm-3\n"
              "RhoBar: %g gcm-3, RhoProBar: %g cm-3\n",
            rho_min,
            rho_min * g2c.g / CUBE( g2c.cm ),
            rho_max,
            rho_max * g2c.g / CUBE( g2c.cm ),
            rho_min * g2c.g / CUBE( g2c.cm ) / cuc.m_p,
            rho_max * g2c.g / CUBE( g2c.cm ) / cuc.m_p,
            RhoBaryon * g2c.g / CUBE( g2c.cm ),
            RhoBaryon * g2c.g / CUBE( g2c.cm ) / cuc.m_p
            );
#ifdef READMACH
    writelog( "[MachNumber], min: %g, max: %g\n", Mmin, Mmax );
#endif
}

