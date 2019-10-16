#include "allvars.h"

void gas_ratio() {
#ifdef GASRATIO
    double m_tot, m_diffuse_cool, m_diffuse_warm, m_diffuse_hot, m_dense, dens, T, m,
           *m_diffuse_cool_all,
           *m_diffuse_warm_all,
           *m_diffuse_hot_all,
           *m_dense_all, *z_all,
           T7, T5, D3;
    long p;
    int i;
    FILE *fd;

    m_tot = m_diffuse_cool = m_diffuse_warm = m_diffuse_hot = m_dense = 0;

    T5 = 1e5;
    T7 = 1e7;
    D3 = 1e3;

    for( p=0; p<N_Gas; p++ ) {
        m = P[p].Mass;
        dens = SphP[p].Density / Time3 / RhoBaryon;
        T = SphP[p].Temp;

        m_tot += m;
        
        if ( T<T5 ) {
            if ( dens>=D3 )
                m_dense += m;
            else
                m_diffuse_cool += m;
            continue;
        }

        if ( T<T7 && T>=T5 ) {
           m_diffuse_warm += m;
            continue;
        }

        if ( T>=T7 ){
            m_diffuse_hot += m;
            continue;
        }
    }

    m_diffuse_cool /= m_tot;
    m_diffuse_warm /= m_tot;
    m_diffuse_hot /= m_tot;
    m_dense /= m_tot;


    if ( ThisTask_Master == 0 ) {
        mymalloc2( z_all, sizeof(double) * NTask );
        mymalloc2( m_diffuse_cool_all, sizeof(double) * NTask );
        mymalloc2( m_diffuse_warm_all, sizeof(double) * NTask );
        mymalloc2( m_diffuse_hot_all, sizeof(double) * NTask );
        mymalloc2( m_dense_all, sizeof(double) * NTask );
    }

    MPI_Gather( &Redshift, 1, MPI_DOUBLE, z_all, 1, MPI_DOUBLE, 0, MpiComm_Master );
    MPI_Gather( &m_diffuse_cool, 1, MPI_DOUBLE, m_diffuse_cool_all, 1, MPI_DOUBLE, 0, MpiComm_Master );
    MPI_Gather( &m_diffuse_warm, 1, MPI_DOUBLE, m_diffuse_warm_all, 1, MPI_DOUBLE, 0, MpiComm_Master );
    MPI_Gather( &m_diffuse_hot, 1, MPI_DOUBLE, m_diffuse_hot_all, 1, MPI_DOUBLE, 0, MpiComm_Master );
    MPI_Gather( &m_dense, 1, MPI_DOUBLE, m_dense_all, 1, MPI_DOUBLE, 0, MpiComm_Master );

    if ( ThisTask_Master == 0 ) {
        fd = myfopen( "w", "%s/gas_ratio.dat", OutputDir );
        for( i=0; i<NTask_Master; i++ )
            fprintf( fd, "%g %g %g %g %g\n",
            z_all[i],
            m_diffuse_cool_all[i],
            m_diffuse_warm_all[i],
            m_diffuse_hot_all[i],
            m_dense_all[i]
            );
        fclose( fd );
    }

    if ( ThisTask_Master == 0 ) {
        myfree( z_all );
        myfree( m_diffuse_cool_all);
        myfree( m_diffuse_warm_all);
        myfree( m_diffuse_hot_all);
        myfree( m_dense_all );
    }
#endif
}
