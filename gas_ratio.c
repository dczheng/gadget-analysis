#include "allvars.h"

void gas_ratio() {
    double m_tot, m_diffuse_cool, m_diffuse_warm, m_diffuse_hot, m_dense, dens, T, m,
           *m_diffuse_cool_all,
           *m_diffuse_warm_all,
           *m_diffuse_hot_all,
           *m_dense_all, *z_all ;
    long p;
    int i;
    char buf[100];
    FILE *fd;

    m_tot = m_diffuse_cool = m_diffuse_warm = m_diffuse_hot = m_dense = 0;

    for( p=0; p<N_Gas; p++ ) {
        m = P[p].Mass;
        dens = SphP[p].Density / All.RhoBaryon;
        T = SphP[p].Temp;

        m_tot += m;
        
        if ( dens >= 1e3 ) {
            m_dense += m;
            continue;
        }

        if ( dens < 1e3 ) {
            if ( T < 1e5 )
                m_diffuse_cool += m;
            if ( T >= 1e5 && T < 1e7 )
                m_diffuse_warm += m;
            if ( T >= 1e7 )
                m_diffuse_hot += m;
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

    MPI_Gather( &All.RedShift, 1, MPI_DOUBLE, z_all, 1, MPI_DOUBLE, 0, MpiComm_Master );
    MPI_Gather( &m_diffuse_cool, 1, MPI_DOUBLE, m_diffuse_cool_all, 1, MPI_DOUBLE, 0, MpiComm_Master );
    MPI_Gather( &m_diffuse_warm, 1, MPI_DOUBLE, m_diffuse_warm_all, 1, MPI_DOUBLE, 0, MpiComm_Master );
    MPI_Gather( &m_diffuse_hot, 1, MPI_DOUBLE, m_diffuse_hot_all, 1, MPI_DOUBLE, 0, MpiComm_Master );
    MPI_Gather( &m_dense, 1, MPI_DOUBLE, m_dense_all, 1, MPI_DOUBLE, 0, MpiComm_Master );

    if ( ThisTask_Master == 0 ) {
        sprintf( buf, "%s/gas_ratio.dat", All.OutputDir );
        fd = fopen( buf, "w" );
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
}
