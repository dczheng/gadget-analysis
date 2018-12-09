#include "allvars.h"

void total_radio_spectrum() {

    long index;
    int Nnu, i, signal, j, k, nu_i;
    double *nu, *flux, numin, numax, dnu, tmp,
           *alist, *fluxlist, *dislist;
    char buf[100];
    FILE *fd;

    if ( ThisTask_Local != 0 )
        return;

    writelog( "total radio spectrum ...\n" );

    Nnu = All.NuNum;
    signal = N_Gas / 10;
    numin = All.NuMin;
    numax = All.NuMax;
    dnu = log( numax/numin) / (Nnu-1);

    mymalloc1( nu, sizeof( double ) * Nnu );
    mymalloc2( flux, sizeof( double ) * Nnu );

    for ( i=0; i<Nnu; i++ ) {
        nu[i] = exp(log(numin) + i * dnu);
    }

    //tmp = 1.0 / ( SQR(All.Time) * ( 4.0 * PI * SQR( All.ComDis * g2c.cm ) ) * ( SQR( All.BoxSize / All.ComDis ) ) );
    tmp = 1.0 / ( ( 4.0 * PI * SQR( All.LumDis * g2c.cm ) ) * ( SQR( All.BoxSize / All.ComDis ) ) );

    if ( All.RedShift > 1e-10 )
            for ( index=0; index<N_Gas; index++ ) {
                    if  ( index % signal == 0 )
                        writelog( "total spectrum [%8li] [%8li] [%5.1f%%]\n",
                               index, N_Gas, ( (double)index )/N_Gas * 100 );

                    for ( i=0; i<Nnu; i++ ) {
                        nu_i = log( nu[i]/All.Time / numin ) / dnu;
                        //tmp = particle_radio( nu[i]*1e6 / All.Time, index );
                        //if ( tmp < 0 )
                         //   endrun( 20181006 );
                         //
                        if ( nu_i < 0 || nu_i >= Nnu )
                            continue;

                        flux[i] += PartRad[index * All.NuNum + nu_i];

                        if ( PartRad[index*All.NuNum+nu_i] * tmp > 10 || PartRad[index*All.NuNum+nu_i] < 0 || flux[i] < 0 ) {
                            printf( "%i, %g, %g, %g\n", nu_i, nu[ nu_i ], PartRad[index*All.NuNum+nu_i], flux[i] );
                            endrun( 20181029 );
                        }

                    }
        }


    for ( i=0; i<Nnu; i++ )
        flux[i] *= tmp;

    create_dir( "TotalSpec" );

    sprintf( buf, "./TotalSpec/Spec_Tot_%.2f.dat", All.RedShift );

    fd = fopen( buf, "w" );

    for( i=0; i<Nnu; i++ )
        fprintf( fd, "%g %g\n", nu[i], flux[i] );

    fclose( fd );

    do_sync_master( "calc total spectrum" );

    //endrun( 20181029 );

    mymalloc1( alist, sizeof( double ) * NTask );
    mymalloc1( dislist, sizeof( double ) * NTask );
    mymalloc1( fluxlist, sizeof( double ) * NTask * Nnu );

    MPI_Gather( &All.Time, 1, MPI_DOUBLE, alist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Gather( flux, Nnu, MPI_DOUBLE, fluxlist, Nnu, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    if ( ThisTask == 0 ) {
        for ( i=0; i<NTask-1; i++ )
            for ( j=i; j<NTask; j++ ) {
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

        for ( i=0; i<NTask-1; i++ )
            dislist[i] = comoving_distance( alist[i] );

        for ( i=0; i<NTask*Nnu; i++ )
            fluxlist[i] /= All.BoxSize; // unit distance

        memset( flux, 0, sizeof(double)*Nnu );

        for ( i=0; i<Nnu; i++ )
            for ( j=0; j<NTask-1; j++) {
                flux[i] += 0.5 * ( fluxlist[j*Nnu+i] + fluxlist[(j+1)*Nnu+i] ) * ( dislist[j+1] - dislist[j] );
            }

        for( i=0; i<Nnu; i++ )
            flux[i] += fluxlist[i] * dislist[0];  // low redshift

        sprintf( buf, "./TotalSpec/Spec_Tot.dat" );
        fd = fopen( buf, "w" );

        for ( i=0; i<Nnu; i++ )
            fprintf( fd, "%g %g\n", nu[i], flux[i] );
        fclose( fd );
    }

    myfree( alist );
    myfree( dislist );
    myfree( fluxlist );
    myfree( nu );
    myfree( flux );

    writelog( "total radio spectrum ... done.\n" );

}
