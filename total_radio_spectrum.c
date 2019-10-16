#include "allvars.h"

void total_radio_spectrum() {
#ifdef TOTSPEC

    long index;
    int Nnu, i, j, k, t;
    double *nu, *flux_local, *flux, numin, numax, dnu, tmp,
           *alist, *fluxlist, *dislist;
    FILE *fd;

    writelog( "total radio spectrum ...\n" );
    Nnu = All.NuNum;
    numin = All.NuMin;
    numax = All.NuMax;
    dnu = log( numax/numin) / (Nnu-1);

    mymalloc1( nu, sizeof( double ) * Nnu );
    mymalloc2( flux_local, sizeof( double ) * Nnu );

    for ( i=0; i<Nnu; i++ ) {
        nu[i] = exp(log(numin) + i * dnu);
    }

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
            //nu_i = log( nu[i]/Time / numin ) / dnu;
            //tmp = particle_radio( nu[i]*1e6 / Time, index );
            //if ( tmp < 0 )
             //   endrun( 20181006 );
             //
             /*
            if ( nu_i < 0 || nu_i >= Nnu )
                continue;
                */

            flux_local[i] += get_particle_radio_index(index, i);

            if ( get_particle_radio_index(index,i) * tmp > 10 ||
                    get_particle_radio_index(index, i) < 0 ||
                    flux_local[i] < 0 ) {
                printf( "%i, %g, %g, %g\n",
                        i, nu[ i ], get_particle_radio_index(index, i), flux_local[i] );
                endrun( 20181029 );
            }

        }
    }

    if ( ThisTask_Local == 0 )
        mymalloc2( flux, sizeof( double ) * Nnu );

    MPI_Reduce( flux_local, flux, Nnu, MPI_DOUBLE, MPI_SUM, 0, MpiComm_Local );
    myfree( nu );
    myfree( flux_local );

    if ( ThisTask_Local )
        return;

    for ( i=0; i<Nnu; i++ )
        flux[i] *= tmp;

    create_dir( "%sTotalSpec", OutputDir );

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
                if ( alist[i] < alist[j] ) {

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

        fd = myfopen( "w", "%s/TotalSpec/Spec_Tot_Comp.dat", OutputDir );
        fprintf( fd, "0 " );
        for( i=0; i<NTask_Master; i++ ) 
            fprintf( fd, "%g ", alist[i] );
        fprintf( fd, "\n" );

        fprintf( fd, "0 " );
        for( i=0; i<NTask_Master; i++ ) 
            fprintf( fd, "%g ", dislist[i] );
        fprintf( fd, "\n" );

        for ( i=0; i<Nnu; i++ ) {
            fprintf( fd, "%g ", numin * exp( i * dnu ) );
            for ( j=0; j<NTask_Master; j++) {
                fprintf( fd, "%g ", fluxlist[j*Nnu+i] );
            }
            fprintf( fd, "\n" );
        }
        fclose( fd );

        for ( i=0; i<NTask_Master*Nnu; i++ ) {
            fluxlist[i] /= BoxSize; // unit distance
        }

        memset( flux, 0, sizeof(double)*Nnu );
        for ( i=0; i<Nnu; i++ )
            for ( j=0; j<NTask_Master-1; j++) {
                flux[i] += 0.5 * ( fluxlist[j*Nnu+i] + fluxlist[(j+1)*Nnu+i] ) * ( dislist[j+1] - dislist[j] );
            }

        for( i=0; i<Nnu; i++ ) {
            flux[i] += fluxlist[i] * dislist[0];  // low redshift
        }

        fd = myfopen( "w", "%s/TotalSpec/Spec_Tot.dat", OutputDir );

        for ( i=0; i<Nnu; i++ )
            fprintf( fd, "%g %g\n", exp(log(numin) + i * dnu), flux[i] );
        fclose( fd );
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
