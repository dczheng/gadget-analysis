#include "allvars.h"

void total_radio_spectrum() {
#ifdef TOTSPEC

    long index;
    int Nnu, i, j, k, t;
    double nu, *flux_local, *flux, numin, numax, dnu, tmp,
           *alist, *fluxlist, *dislist;
    FILE *fd;

    writelog( "total radio spectrum ...\n" );
    Nnu = All.FreqN;
    numin = All.FreqMin;
    numax = All.FreqMax;
    dnu = log( numax/numin) / (Nnu-1);

    mymalloc2( flux_local, sizeof( double ) * Nnu );
    tmp = 1.0 / ( ( 4.0 * PI * SQR( LumDis * g2c.cm ) ) * BoxSolidAngle ) ;

    t = N_Gas / 10;
    for ( index=0; index<N_Gas; index++ ) {
        if  ( index % t == 0 ) {
                writelog( "total spectrum [%8li] [%8li] [%5.1f%%]\n",
                       index, N_Gas, ( (double)index )/N_Gas * 100 );
        }

        if ( index % NTask_Local != ThisTask_Local )
            continue;

        for ( i=0; i<Nnu; i++ ) {

            nu = exp(log(numin) + i * dnu);
            flux_local[i] += get_particle_radio(index, nu);

            if ( get_particle_radio(index,nu) * tmp > 10 ||
                    get_particle_radio(index, nu) < 0 ||
                    flux_local[i] < 0 ) {
                printf( "%i, %g, %g, %g\n",
                        i, nu, get_particle_radio(index, nu),
                        flux_local[i] );
                endrun( 20181029 );
            }
        }
    }

    if ( ThisTask_Local == 0 )
        mymalloc2( flux, sizeof( double ) * Nnu );

    MPI_Reduce( flux_local, flux, Nnu, MPI_DOUBLE, MPI_SUM, 0, MpiComm_Local );
    myfree( flux_local );

    create_dir( "%sTotalSpec", OutputDir );

    if ( ThisTask_Local )
        return;

    for ( i=0; i<Nnu; i++ )
        flux[i] *= tmp;

    if ( ThisTask_Master == 0 ) {
        mymalloc2( alist, sizeof( double ) * NTask_Master );
        mymalloc2( dislist, sizeof( double ) * NTask_Master );
        mymalloc2( fluxlist, sizeof( double ) * NTask_Master * Nnu );
    }

    MPI_Gather( &Time, 1, MPI_DOUBLE, alist, 1, MPI_DOUBLE, 0, MpiComm_Master );
    MPI_Gather( flux, Nnu, MPI_DOUBLE, fluxlist, Nnu, MPI_DOUBLE, 0, MpiComm_Master );

    if ( ThisTask_Master == 0 ) {
        for ( i=0; i<NTask_Master-1; i++ )
            for ( j=i; j<NTask_Master; j++ ) {
                if ( alist[i] > alist[j] ) {

                        tmp = alist[i];
                        alist[i] = alist[j];
                        alist[j] = tmp;

                        for ( k=0; k<Nnu; k++ ) {
                            tmp = fluxlist[ i*Nnu + k ];
                            fluxlist[ i*Nnu + k ] = fluxlist[ j*Nnu + k ];
                            fluxlist[ j*Nnu + k ] = tmp;
                        }
                }
            }

        for ( i=0; i<NTask_Master; i++ ) {
            dislist[i] = comoving_distance( alist[i] );
        }

#define myprint1( _fd ) {\
        int _i;\
        fprintf( _fd, "0 " );\
        for( _i=0; _i<Nnu; _i++ ) \
            fprintf( _fd, "%g ", numin * exp(_i*dnu)*Time/1e6 );\
        fprintf( _fd, "\n" );\
}

#define myprint2( _fd, _d, _idx ) {\
        int _i;\
        fprintf( _fd, "%g ", 1/alist[_idx]-1 );\
        for( _i=0; _i<Nnu; _i++ ) \
            fprintf( _fd, "%g ", _d[_i] );\
        fprintf( _fd, "\n" );\
}

        fd = myfopen( "w", "%s/TotalSpec/Spec_Tot_Comp.dat", OutputDir );
        myprint1( fd );
        for( i=0; i<NTask_Master; i++ ) {
            myprint2( fd, (fluxlist+i*Nnu), i );
        }
        fclose( fd );

        for ( i=0; i<NTask_Master*Nnu; i++ ) {
            fluxlist[i] /= BoxSize; // unit distance
        }

        fd = myfopen( "w", "%s/TotalSpec/Spec_Tot.dat", OutputDir );
        myprint1( fd );

        for( k=1; k<NTask_Master; k++ ) {
            memset( flux, 0, sizeof(double)*Nnu );
            for ( i=0; i<Nnu; i++ )
                for ( j=1; j<=k; j++) {
                    flux[i] += 0.5 * ( fluxlist[(j-1)*Nnu+i] + fluxlist[j*Nnu+i] )
                                   * ( dislist[j-1] - dislist[j] );
                }
            myprint2( fd, flux, k );
        }

        fclose( fd );
#undef myprint1
#undef myprint2
    }

    if ( ThisTask_Master == 0 ) {
        myfree( alist );
        myfree( dislist );
        myfree( fluxlist );
    }

    myfree( flux );
    writelog( "total radio spectrum ... done.\n" );

#endif
}
