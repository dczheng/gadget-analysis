#include "allvars.h"
#include "drfftw.h"
#ifdef GROUP

int group_present( long index ) {

    if ( Gprops[index].mass >= All.GroupMassMin) 
            return 1;
    else {

        return 0;
    }

}

inline double get_group_size( struct group_properties *g ) {
    return vmax( g->size[proj_i], g->size[proj_j] );
}

#ifdef GROUPPOT
double group_pot( long index ) {

    long p;
    double e;
    struct group_properties *g;

    g = &Gprops[index];
    p = g->Head;
    e = 0;
    while( p >= 0 ) {
        e += P[p].Pot * P[p].Mass;
        p = FoFNext[p];
    }
    return e;

}
#endif
//#define GROUP_POT2_TEST
#ifdef GROUP_POT2_TEST
double group_pot2_direct( long index ) {
    /*
    for test
    */

    double *m, *pos, r, e;
    long p;
    int i, j, k, N;
    struct group_properties *g;
    g = &Gprops[index];

    printf( "npart: " );
    for( i=0; i<6; i++ )
        printf( "%li ", g->npart[i] );
    printf( "\n" );
    N =  g->npart[0];
    printf( "N: %i\n", N );
    mymalloc2( m, N * sizeof(double) );
    mymalloc2( pos, N * sizeof(double) * 3 );

    index = 0;
    p = g->Head;
    while( p >= 0 ) {
        if ( P[p].Type != 0 ) {
            p = FoFNext[p];
            continue;
        }

        m[index] = P[p].Mass;

        for( i=0; i<3; i++ )
            pos[index*3+i] = P[p].Pos[i];
        index ++;
        p = FoFNext[p];
    }

    e = 0;
    for( i=0; i<N-1; i++ ) {
        for( j=i+1; j<N; j++ ) {
            r = 0;
            for(k=0; k<3; k++)
                r += SQR( PERIODIC_HALF( pos[i*3+k] - pos[j*3+k] ) );
            r = sqrt( r );
            e += -G * m[i] * m[j] / r;
        }
        printf( "N: %i, i: %i, e: %g\n", N, i, e );
    }

    myfree( m );
    myfree( pos );

    return e;
}
#endif

#ifdef GROUPPOT
double group_pot2( long index ) {

    long p, pp;
    int i, j, k, m, NGrid, NGrid3, ks[3], k2, ip, ngbnum, tabindex;
    double e, L, *data, fac, f[3], smth, ff, asmth2, asmth, rcut, r, asmthfac;
    struct group_properties *g;
#ifdef GROUP_POT2_TEST
    index = 2;
#endif

    fftw_real *rhogrid, *massgrid;
    fftw_complex *fft_of_rhogrid;
    rfftwnd_plan fft_plan, fft_plan_inv;
    NGrid = All.GroupPotGrid;
    NGrid3 = 2 * ( NGrid / 2 + 1  );

    mymalloc2( rhogrid, SQR(NGrid) * NGrid3 * sizeof(fftw_real)  );
    mymalloc2( massgrid, SQR(NGrid) * NGrid3 * sizeof(fftw_real)  );
    fft_of_rhogrid = ( fftw_complex*  )rhogrid;

    g = &Gprops[index];
    L = g->size[0];
    for( i=1; i<3; i++ )
        vmax2( L, g->size[i] );
    L *= 1.01;
    L *= 2;

    k = 0;
    for( i=0; i<6; i++ )
        if ( g->npart[i]>0 ) 
            if ( ((All.TreePartType>>i)&1) == 0 )
                endruns( "error" );

    mymalloc1( Ngblist, g->Len * sizeof(long) );

    asmth = ASMTH * ( L/NGrid );
    rcut = RCUT * asmth;
    asmth2 = ( 2*PI ) / L * asmth;
    asmth2 *= asmth2;
    asmthfac = 0.5 / asmth * ( NSRPTAB/3.0 );

    fft_plan = rfftw3d_create_plan( NGrid, NGrid, NGrid,
            FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE ); 
    fft_plan_inv = rfftw3d_create_plan( NGrid, NGrid, NGrid,
            FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE ); 
    mymalloc2( data, sizeof( double ) * g->Len * 4  );
    p = g->Head;
    ip = 0;
#ifdef GROUP_POT2_TEST
    double mtot=0, mtot2=0;
#endif
    while( p >= 0 ) {
#ifdef GROUP_POT2_TEST
        if ( P[p].Type != 0 ) {
            p = FoFNext[p];
            continue;
        }
#endif
        for( i=0; i<3; i++ )
            data[ip*4+i] = PERIODIC_HALF( P[p].Pos[i] - g->cm[i] ) + L/2;

#ifdef GROUP_POT2_TEST
       mtot += P[p].Mass;
#endif
        data[ip*4+3] = P[p].Mass;
        ip ++;
        p = FoFNext[p];
    }

    field_to_grid_fft( data, rhogrid, L, NGrid, ip, 1 );

#ifdef GROUP_POT2_TEST
    for( i=0; i<NGrid; i++ )
        for( j=0; j<NGrid; j++ )
            for( k=0; k<NGrid; k++ ) {
                ip = i * NGrid * NGrid3 + j * NGrid3 + k;
                mtot2 += rhogrid[ip];
            }
    printf( "mtot: %g, mtot2: %g\n", mtot, mtot2 );
#endif

    memcpy( massgrid, rhogrid, sizeof(fftw_real) * SQR(NGrid) * NGrid3 );
    myfree( data );

    fac = G / ( PI * L );
    /*
     to get potential
     [4 * PI * G / ( k_i * 2*Pi / L )^2]  * [1 / ( L/N )^3] * N^3,
     1 / ( L/N )^3: mass to density;
     N^3: implementation of fftw, that is FFTW computes an unnormalized transform,
    */

#ifdef GROUP_POT2_TEST
    for( i=0; i<100; i++ )
        printf( "%g ", rhogrid[NGrid/2 * NGrid * NGrid3 + NGrid/2 * NGrid3 + i] );
    printf( "\n" );
#endif

    rfftwnd_one_real_to_complex( fft_plan, rhogrid, NULL  );

    for( i=0; i<NGrid; i++ )
        for( j=0; j<NGrid; j++ )
            for( k=0; k<NGrid/2+1; k++ ) {

                if ( i==0 && j==0 && k==0 ) {
                    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0;
                    continue;
                }

                ks[0] = ( i>NGrid/2 ) ? i-NGrid : i;
                ks[1] = ( j>NGrid/2 ) ? j-NGrid : j;
                ks[2] = ( k>NGrid/2 ) ? k-NGrid : k;

                for( m=0, k2=0; m<3; m++ ) {
                    k2 += ks[m]*ks[m];
                    f[m] = 1;
                    if ( ks[m] != 0 ) {
                        f[m] = ( PI * ks[m] ) / NGrid;
                        f[m] = sin(f[m]) / f[m];
                    }
                }

                ff = 1 / ( f[0] * f[1] * f[2] );
                smth = - exp( -k2 * asmth2 ) / k2 * fac * ff * ff * ff * ff;

                ip = i * NGrid * ( NGrid/2+1 ) + j * ( NGrid/2+1 ) + k;
                fft_of_rhogrid[ ip ].re *= smth; 
                fft_of_rhogrid[ ip ].im *= smth;
            }

    rfftwnd_one_complex_to_real( fft_plan_inv, fft_of_rhogrid, NULL  );

#ifdef GROUP_POT2_TEST
    for( i=0; i<100; i++ )
        printf( "%g ", rhogrid[NGrid/2 * NGrid * NGrid3 + NGrid/2 * NGrid3 + i] );
    printf( "\n" );
#endif

    e = 0;
    for( i=0; i<NGrid; i++ )
        for( j=0; j<NGrid; j++ )
            for( k=0; k<NGrid; k++ ) {
                ip = i * NGrid * NGrid3 + j * NGrid3 + k;
                e += rhogrid[ip] * massgrid[ip];
            }

/*
    p = g->Head;
    ip = 0;
    e = 0;
    fac = NGrid / L;
    while( p >= 0 ) {

#ifdef GROUP_POT2_TEST
        if ( P[p].Type != 0 ) {
            p = FoFNext[p];
            continue;
        }
#endif
         ip = 
            (PERIODIC_HALF( P[p].Pos[0] - g->cm[0] ) + L/2) * fac * NGrid * NGrid3 +
            (PERIODIC_HALF( P[p].Pos[1] - g->cm[1] ) + L/2) * fac * NGrid3 +
            (PERIODIC_HALF( P[p].Pos[2] - g->cm[2] ) + L/2) * fac ;
         e += P[p].Mass * rhogrid[ip];
        p = FoFNext[p];
    }
*/

/* now, compute shortrange potential energy */
    for( p=0; p<NumPart; p++ )
        P[p].flag = 0;

    p = g->Head;
    while( p >= 0 ) {
#ifdef GROUP_POT2_TEST
        if ( P[p].Type != 0 ) {
            p = FoFNext[p];
            continue;
        }
#endif
        P[p].flag = 1;
        p = FoFNext[p];
    }

    p = g->Head;
    while( p >= 0 ) {
#ifdef GROUP_POT2_TEST
        if ( P[p].Type != 0 ) {
            p = FoFNext[p];
            continue;
        }
#endif
        ngbnum = ngb( P[p].Pos, rcut, 1 );

        //printf( "%i\n", ngbnum );
        for( pp=0; pp<ngbnum; pp++ ) {
            r = 0;
            for(k=0; k<3; k++)
                r += SQR( PERIODIC_HALF( P[p].Pos[k] - P[pp].Pos[k] ) );
            r = sqrt(r);
            tabindex = (int)( asmthfac * r );
            //printf( "%i %i\n", tabindex, NSRPTAB );
            if ( tabindex < NSRPTAB )
                e += -0.5 * G * P[p].Mass * P[pp].Mass / r *
                ShortRangeTablePotential[tabindex];
        }

        p = FoFNext[p];
    }


#ifdef GROUP_POT2_TEST
    FILE *fd;
    double e_direct, t;
    fd = myfopen( "w", "group-pot2-test-pot.dat" );
    for( k=0; k<NGrid; k++ ) {
        for( j=0; j<NGrid; j++ ) {
            t = 0;
            for( i=0; i<NGrid; i++ ) {
                ip = i * NGrid * NGrid3 + j * NGrid3 + k;
                t += rhogrid[ip];
            }
            fprintf( fd, "%g ", t );
        }
        fprintf( fd, "\n" );
   }
   fclose( fd );

    fd = myfopen( "w", "group-pot2-test-mass.dat" );
    for( k=0; k<NGrid; k++ ) {
        for( j=0; j<NGrid; j++ ) {
            t = 0;
            for( i=0; i<NGrid; i++ ) {
                ip = i * NGrid * NGrid3 + j * NGrid3 + k;
                t += massgrid[ip];
            }
            fprintf( fd, "%g ", t );
        }
        fprintf( fd, "\n" );
   }
   fclose( fd );
   //e_direct = group_pot2_direct( index );
   printf( "pot fft: %g \n", e );
   printf( "pot direct: %g \n", e_direct );

   endruns( "group-pot2-test" );
#endif

    rfftwnd_destroy_plan( fft_plan  );
    rfftwnd_destroy_plan( fft_plan_inv );
    myfree( rhogrid );
    myfree( massgrid );
    myfree( Ngblist );

    return e;

}
#endif

#ifdef GROUPKIN
double group_kin( long index ) {

    long p, i;
    double e, m;
    struct group_properties *g;

    g = &Gprops[index];

    p = g->Head;
    e = 0;
    while( p >= 0 ) {
        e += 0.5 * P[p].Mass * 
        ( SQR( P[p].Vel[0] ) + 
                    SQR( P[p].Vel[1] ) + 
                    SQR( P[p].Vel[2] ) ) * Time;
        p = FoFNext[p];

    }
    for( i=0,m=0; i<6; i++ )
        m += g->mass_table[i];
    e -= 0.5 * m *  ( 
                SQR( g->vel[0] ) +
                SQR( g->vel[1] ) +
                SQR( g->vel[2] )
            ) * Time;

    return e;

}
#endif

#ifdef GROUPVELDISP
double group_vel_disp( long index ) {

    long p;
    double v_disp, v_mean, t;
    struct group_properties *g;

    g = &Gprops[index];

    p = g->Head;
    v_mean = 0;
    while( p >= 0 ) {

        if ( P[p].Type == 0 ) {
            t = sqrt( SQR( P[p].Vel[0] ) + 
                        SQR( P[p].Vel[1] ) + 
                        SQR( P[p].Vel[2] ) );
            v_mean += t * sqrt( Time );
        }
        p = FoFNext[p];

    }
    v_mean /= g->npart[0];
    //printf( "v_mean: %g\n", v_mean );

    p = g->Head;
    v_disp = 0;
    while( p >= 0 ) {

        if ( P[p].Type == 0 ) {
            t = sqrt( SQR( P[p].Vel[0] ) + 
                        SQR( P[p].Vel[1] ) + 
                        SQR( P[p].Vel[2] ) );
            t =  t * sqrt( Time ) - v_mean;
            v_disp += t*t;
        }
        p = FoFNext[p];

    }
    v_disp /= g->npart[0];
    v_disp = sqrt( v_disp );

    return v_disp;

}
#endif

#ifdef GROUPVIRIAL
void group_virial() {
     FILE *fd;
     int index; 
     double v_disp, ek, ep, ep2;
     put_header( "group virial" );
     create_dir( "%sVirial", GroupDir  );
     if ( ThisTask_Local == 0 )
        fd = myfopen( "w", "%sVirial/Virial_%03i.dat",
            GroupDir, SnapIndex);
     for ( index=0; index<Ngroups; index++  ) {
         if ( !group_present( index  )  )
             break;
             writelog( "group: %i\n", index );
             v_disp = group_vel_disp( index );
             ep = group_pot( index );
             ek = group_kin( index );
             ep2 = group_pot2( index );
             if ( ThisTask_Local == 0 )
                fprintf( fd, "%g %g %g %g\n", ep, ek, v_disp, ep2 );
     }
     if ( ThisTask_Local == 0 )
        fclose( fd );
     put_end();
}
#endif


/*
void output_group( int gidx ) {
    long p, i;
    FILE *fd;
    fd = myfopen( "w", "g_%04i.dat", gidx );
    p = Gprops[gidx].Head;
    while( p>=0 ) {
        if ( P[p].Type != 0 ) {
            p = FoFNext[p];
            continue;
        }

        if ( SphP[p].CRE_C == 0 ) {
            p = FoFNext[p];
            continue;
        }

        fprintf( fd, "%e %5.2f %e %e %e ",
        SphP[p].CRE_C,
        SphP[p].CRE_Alpha,
        SphP[p].CRE_qmin,
        SphP[p].CRE_qmax,
        get_B( p ) * 1e6
        );
        for( i=0; i<All.NuNum; i++ ) {
            fprintf( fd, "%e ", get_particle_radio_index( p, i ) );
        }
        fprintf( fd, "\n" );
        p = FoFNext[p];
    }
    fclose( fd );
    //hendruns( "output_group" );

}
*/

void group_analysis() {

    int proj_tmp[3], i;
    char Sproj_tmp;

    put_header( "group analysis" );

    //output_group( 0 );

#ifdef GROUPTOTLUM
    group_tot_lum();
#endif

#ifdef GROUPVIRIAL
    group_virial();
#endif

#ifdef GROUPSPEC
        group_spectrum();
        //group_spectrum_index();
#endif

#ifdef GROUPELECSPEC
        group_electron_spectrum();
        put_sep0;
#endif

    for( i=0; i<3; i++ )
        proj_tmp[i] = Proj[i];
    Sproj_tmp = Sproj;

    for( i=0; i<3; i++ ) {
        proj_k = i;
        proj_i = ( i + 1  ) % 3;
        proj_j = ( i + 2  ) % 3;
        Sproj = i + 'x';
        group_plot();
    }

    for( i=0; i<3; i++ )
        Proj[i] = proj_tmp[i];
    Sproj = Sproj_tmp;


#ifdef GROUPTEMPPROFILE
        group_temp_profile();
#endif

#ifdef GROUPTEMPSTACK
        group_temp_stack();
#endif

#ifdef GROUPGASRATIO
        group_gas_ratio();
#endif

    reset_img();
    put_sep0;
    put_end();

}

#endif


