#include "allvars.h"

void corr_Tdiff_dens() {

    double z[2], *T, *Dens, dx, *T0, *Dens0, Tmin, Tmax, Densmin, Densmax,
          // Tmm[4], Densmm[4],
           *Tdiff, *Densdiff, Tdiffbar, Densbar,
           *corr_local, r, rmin, rmax, dr, rr, *corr;
    int N, N2, N3, i, j, k, index, x[3], *num, ri, index0, l,
        i0, j0, k0, rN;
    MPI_Status status;
    long p;
    FILE *fd;

    writelog( "Correlation of T difference and density ...\n" );

    if ( ThisTask_Local == 0 ) {
        MPI_Gather( &All.RedShift, 1, MPI_DOUBLE, z, 1, MPI_DOUBLE, 0, MpiComm_Master );
        if ( ThisTask_Master == 0 ) {
            if ( fabs( z[0] - z[1] ) > 1e-5 ) {
                printf( "The same redshift is required by corralation of "
                    "T difference and density\n" );
                endrun( 20190720 );
            }
        }
    }

    N = All.CorrGrid;
    N2 = SQR(N);
    N3 = N2*N;
    dx = All.BoxSize / N;

    rN = All.CorrRN;
    rmin = All.CorrRMin;
    rmax = All.CorrRMax;

    if ( rmin < dx ) {
        writelog( "rmin (%g) < dx (%g)!\n", rmin, dx );
        while( rmin < dx )
            rmin *= 1.001;
        writelog ( "set rmin to: %g\n", rmin );
    }
    if ( rmax < dx ) {
        writelog( "rmax (%g) < dx (%g)!\n", rmax, dx );
        endrun( 20190722 );
    }


    dr = ( rmax - rmin ) / ( rN - 1 );

    if ( dr < rmin ) {
        writelog( "reset `dr` and `rN`\n" );
        dr = dx;
        rN = ( rmax - rmin ) / dr + 1;
    }


    if ( ThisTask_Local == 0 ) {
        mymalloc2( T, sizeof( double ) * N3 );
        mymalloc2( num, sizeof( int ) * N3 );

        if ( ThisTask_Master == 0 ) {
            mymalloc2( T0, sizeof( double ) * N3 );
            mymalloc2( Dens0, sizeof( double ) * N3 );
            mymalloc2( Densdiff, sizeof( double ) * N3 );
        }
    }

    mymalloc2( Dens, sizeof( double ) * N3 );
    mymalloc2( Tdiff, sizeof( double ) * N3 );
    mymalloc2( corr_local, sizeof( double ) * rN );


    if ( ThisTask_Local == 0 ) {

        Tmin = Densmin = DBL_MAX;
        Tmax= Densmax = -DBL_MAX;

        for ( p=0; p<N_Gas; p++ ) {
            Tmin = ( SphP[p].Temp < Tmin ) ? SphP[p].Temp : Tmin;
            Tmax = ( SphP[p].Temp > Tmax ) ? SphP[p].Temp : Tmax;

            Densmin = ( SphP[p].Density < Densmin ) ? SphP[p].Density : Densmin;
            Densmax = ( SphP[p].Density < Densmax ) ? SphP[p].Density : Densmax;

            for( i=0; i<3; i++ ) {
                x[i] = P[p].Pos[i] / dx;
                x[i] = x[i] % N;
            }
            index = x[0]*N*N + x[1]*N + x[2];
            T[index] += SphP[p].Temp;
            Dens[index] += SphP[p].Density;
            num[index] ++;
        }

        /*
        for( index=0; index<N3; index++ ) {
            if ( num[index] > 0 ) {
                T[index] /= num[index];
                Dens[index] /= num[index];
            }
        }
        */

        /*
        MPI_Gather( &Tmin, 1, MPI_DOUBLE, Tmm, 1, MPI_DOUBLE, 0, MpiComm_Master );
        MPI_Gather( &Tmax, 1, MPI_DOUBLE, Tmm+2, 1, MPI_DOUBLE, 0, MpiComm_Master );

        t = Tmm[1];
        Tmm[1] = Tmm[2];
        Tmm[2] = t;

        MPI_Gather( &Densmin, 1, MPI_DOUBLE, Densmm, 1, MPI_DOUBLE, 0, MpiComm_Master );
        MPI_Gather( &Densmax, 1, MPI_DOUBLE, Densmm+2, 1, MPI_DOUBLE, 0, MpiComm_Master );

        t = Densmm[1];
        Densmm[1] = Densmm[2];
        Densmm[2] = t;
        */

    }


    if ( ThisTask_Local == 0 ) {

        if ( ThisTask_Master != 0 ) {
            MPI_Send( T, N3, MPI_DOUBLE, 0, 0, MpiComm_Master );
            MPI_Send( Dens, N3, MPI_DOUBLE, 0, 0, MpiComm_Master );
        }
        else {
            MPI_Recv( T0, N3, MPI_DOUBLE, 1, 0, MpiComm_Master, &status);
            MPI_Recv( Dens0, N3, MPI_DOUBLE, 1, 0, MpiComm_Master, &status);

            Tdiffbar = Densbar = 0;
            for( index=0; index<N3; index++ ) {

                if ( T0[index] > 0 ) {
                    Tdiff[index] = ( T[index] - T0[index] ) / T0[index];
                }

                if ( Dens0[index] > 0 ) {
                    Densdiff[index] = ( Dens[index]-Dens0[index] ) / Dens0[index];
                }

                Tdiffbar += Tdiff[index];
                Densbar += Dens[index];

            }

            Tdiffbar /= N3;
            Densbar /= N3;

            FILE  *fd;
            fd = fopen( "T_Dens.dat", "w" );
            for( i=0; i<N; i++ )
                for( j=0; j<N; j++ )
                    for( k=0; k<N; k++ ) {
                        index = i * N2 + j * N + k;
                        fprintf( fd, "%g %g %g %g %g %g %g %g %g\n",
                            i * dx, j * dx, k * dx,
                            T[index], T0[index], Tdiff[index],
                            log10(Dens[index]), log10(Dens0[index]),
                            Densdiff[index]
                        );
            }

            fclose( fd );

            for( index=0; index<N3; index++ ) {
                Tdiff[index] = ( Tdiff[index] - Tdiffbar ) / Tdiffbar;
                Dens[index] = ( Dens[index] - Densbar ) / Densbar;
            }

            /*
            fd = fopen( "Dens.dat", "w" );
            for( p=0; p<N_Gas; p++ ) {
                fprintf( fd, "%g %g %g %g\n",
                        P[p].Pos[0],
                        P[p].Pos[1],
                        P[p].Pos[2],
                        log10(SphP[p].Density)
                        );
            }
            fclose( fd );
            */
        }
    }

    MPI_Bcast( Tdiff, sizeof(double)*N3, MPI_BYTE, 0, MPI_COMM_WORLD );
    MPI_Bcast( Dens, sizeof(double)*N3, MPI_BYTE, 0, MPI_COMM_WORLD );

    //printf( "%i\n", ThisTask );

    for( ri=1; ri<rN; ri++ ) {
        r = rmin + ri * dr;
        writelog( "%i, %g\n", ri, r );
        l = r / dx + 1;

        for( i=0; i<N; i++ ) {
            for( j=0; j<N; j++ ) {
                for ( k=0; k<N; k++ ) {
                    index = i * N2 + j * N + k;

                    if ( index % NTask != ThisTask )
                        continue;

                    for( i0=i-l; i0<=i+l; i0++ )
                        for( j0=j-l; j0<=j+l; j0++ )
                            for( k0=k-l; k0<=k+l; k0++ ) {

                                if ( i0 < 0 ||
                                     i0 > N-1  ||
                                     j0 < 0  ||
                                     j0 > N-1  ||
                                     k0 < 0  ||
                                     k0 > N-1 )
                                    continue;

                                rr = sqrt(
                                        SQR( i0-i ) +
                                        SQR( j0-j ) +
                                        SQR( k0-k )
                                        ) * dx;

                                if ( fabs( r-rr ) > dx )
                                    continue;

                                index0 = i0*N2 + j0*N + k0;

                                corr_local[ri] += Tdiff[index] * Dens[index0]
                                    * CUBE( dx );
                                //corr_local[ri] += Tdiff[index] * Dens[index0];
                            }
                }
            }
        }
    }

    if ( ThisTask == 0 )
        mymalloc2( corr, sizeof( double ) * rN );

    MPI_Reduce( corr_local, corr, rN, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    if ( ThisTask == 0 ) {
        fd = fopen( "corr_Tdiff_dens.dat", "w" );
        for( ri=0; ri<rN; ri++ ) {
            fprintf( fd, "%g %g\n",
                    rmin + ri * dr, corr[ri]
                    );
        }

        fclose( fd );

        myfree( corr );
    }

    if ( ThisTask_Local == 0 ) {
        myfree( T );
        myfree( num );

        if ( ThisTask_Master == 0 ) {
            myfree( T0 );
            myfree( Dens0 );
            myfree( Densdiff );
        }
    }

    myfree( Dens );
    myfree( Tdiff );
    myfree( corr_local );


}
