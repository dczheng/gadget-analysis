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

#ifdef GROUP_POT_TEST
double group_pot_direct( long index ) {
    /*
    for test
    */

    double *m, *pos, r, e, e_tot;
    long p;
    int i, j, k, N;
    struct group_properties *g;
    g = &Gprops[index];

    N =  g->npart[0];
    if ( ThisTask_Local == 0 ) {
        printf( "npart: " );
        for( i=0; i<6; i++ )
            printf( "%li ", g->npart[i] );
        printf( "\n" );
        printf( "N: %i\n", N );
    }

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
        if ( i % NTask_Local != ThisTask_Local )
            continue;

        for( j=i+1; j<N; j++ ) {
            r = 0;
            for(k=0; k<3; k++)
                r += SQR( PERIODIC_HALF( pos[i*3+k] - pos[j*3+k] ) );
            r = sqrt( r );
            e += -G * m[i] * m[j] / r;
        }
#ifdef GROUP_POT_TEST2
        if ( i%(N/10) == 0 && ThisTask_Local == 0 )
            printf( "%5.2f%%, e: %g\n", ((double)i)/N*100, e );
#endif
    }

    myfree( m );
    myfree( pos );

    MPI_Allreduce( &e, &e_tot, 1, MPI_DOUBLE, MPI_SUM, MpiComm_Local );

    return e_tot;
}
#endif

#ifdef GROUPPOT
//#define TREEPOT
double group_pot( long index ) {

    long p;
    int i, j, k, NGrid, ip, N;
    double e, L, *data, *pot;
    struct group_properties *g;

    NGrid = All.GroupPotGrid;
    fftw_real *massgrid;

    mymalloc2( massgrid, CUBE(NGrid) * sizeof(fftw_real)  );
    mymalloc2( pot, CUBE(NGrid) * sizeof(double)  );

    g = &Gprops[index];
    L = g->size[0];
    for( i=1; i<3; i++ )
        vmax2( L, g->size[i] );
    L *= 2;
    L += 4 * L/NGrid;   
#ifdef GROUP_POT_TEST
    N = g->npart[0];
#else
    N = g->Len;
#endif

    check_fof( 0, 10 );

    mymalloc2( data, sizeof(double)*N*4 );

    for( i=0; i<6; i++ )
        if ( g->npart[i]>0 ) 
            if ( ((All.TreePartType>>i)&1) == 0 )
                endruns( "error" );

    ip = 0;
    p = g->Head;
    while( p >= 0 ) {
        printf( "%i,", P[p].Type );
        if ( P[p].Pos[0]<70000 ) {
            printf( "xx: %g, %li\n",P[p].Pos[0], index );
            endruns( "xxx" );
        }
        for( i=0; i<3; i++ ) {
            printf( "%g ", P[p].Pos[i] );
        }
        for( i=0; i<3; i++ )
            printf( "%g ", g->cm[i] );
        printf( "\n" );
        p = FoFNext[p];
        continue;
#ifdef GROUP_POT_TEST
        if ( P[p].Type ) {
            p = FoFNext[p];
            continue;
        }
#endif
        for( i=0; i<3; i++ ) {
            data[ip*4+i] = PERIODIC_HALF( P[p].Pos[i] - g->cm[i] ) + L/2;
            if ( fabs(PERIODIC_HALF( P[p].Pos[i] - g->cm[i] )) > L/2){
                printf( "%li, %i, %i, %g %g %g\n",
                        index,
                    P[p].Type, i, P[p].Pos[i], g->cm[i], L );
                endruns( "can't occur!" );
            }
        }

        data[ip*4+3] = P[p].Mass;
        ip ++;
        p = FoFNext[p];
    }

    if ( ip != N )
           endruns( "can't occur!" );

    field_to_grid_fft( data, massgrid, L, NGrid, N, 0, 0 );
    pm_potential( data, N, L, NGrid, 0, pot );

/*
    for( i=0; i<10; i++ )
        printf( "%g\n",pot[ NGrid/2*NGrid*NGrid + NGrid/2*NGrid + (NGrid/2-5)+i ]);
*/

    e = 0;
    for( i=0; i<NGrid; i++ )
        for( j=0; j<NGrid; j++ )
            for( k=0; k<NGrid; k++ ) {
                ip = i * NGrid * NGrid + j * NGrid + k;
                e += 0.5 * pot[ip] * massgrid[ip];
            }

/*
    ip = 0;
    e = 0;
    double fac = NGrid / L;
    for( i=0; i<N; i++ ){
            ip =
            data[i*4+0] * fac * NGrid * NGrid +
            data[i*4+1] * fac * NGrid +
            data[i*4+2] * fac;

        e += data[i*4+3] * pot[ip];
    }
*/

#ifdef TREEPOT
/* now, compute shortrange potential energy */
    int tabindex, ngbnum;
    long pp;
    double rcut, asmthfac, tabindex;
    mymalloc1( Ngblist, g->Len * sizeof(long) );

    rcut = asmth;
    //rcut = RCUT * asmth;
    asmthfac = 0.5 / asmth * ( NSRPTAB/3.0 );

    for( p=0; p<NumPart; p++ )
        P[p].flag = 0;

    p = g->Head;
    while( p >= 0 ) {
#ifdef GROUP_POT_TEST
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
#ifdef GROUP_POT_TEST
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
    myfree( Ngblist );
#endif


    myfree( massgrid );
    myfree( pot );
    myfree( data );

    return e;

}
#endif

#ifdef OUTPUTGROUP
#define OUTPUTGROUP_DEBUG
#ifdef OUTPUTGROUP_DEBUG
double group_ek( long index ) {
    long p;
    int i;
    struct group_properties *g;
    double t, vc[3],v2, m, ek, mtot;
    g = &Gprops[index];

    vc[0] = vc[1] = vc[2] = 0;
    ek = 0;
    p = g->Head; 
    while( p>=0 ) {
        m = P[p].Mass;
        for( i=0,v2=0; i<3; i++ ) {
            t = P[p].Vel[i] / sqrt(Time);
            vc[i] += t * m;
            v2 += t*t;
        }
        ek += 0.5 * m * v2;
        p = FoFNext[p];
    }

    for(i=0,mtot=0; i<6; i++)
        mtot += g->mass_table[i];
    for(i=0,v2=0; i<3; i++) {
        vc[i] /= mtot;
        v2 += vc[i]*vc[i];
    }

    ek -= 0.5 * mtot * v2;
    return ek;

}
#endif

void output_group() {

    int index, i;
    long p;
    double ep, ek, m_diffuse, m_warmhot, m_hot, m_condensed, m, d, t;
    FILE *fd;
    struct group_properties *g;
    put_header( "output_group" );
#ifdef OUTPUTGROUPLUM
    double lum;
#endif


    fd = myfopen( "w", "%s/group_%03i.csv",
            GroupDir, SnapIndex);
    fprintf( fd, "mtot,"
            "x,y,z,"
            "mgas,mdm,mstar,mdiffuse,mwarmhot,mhot,mcondensed,"
#ifdef OUTPUTGROUPLUM
           "Lum,"
#endif
           "r200,ek,ep,vr,v_mean,v_disp"
           "\n");

    for ( index=0; index<Ngroups; index++  ) {

        if ( !group_present(index) )
            continue;

        //printf( "%i\n", index );
        g = &Gprops[index];
        p = g->Head; 
        m_diffuse = m_warmhot = m_hot = m_condensed = 0;
        while( p>=0 ) {
            if ( P[p].Type ) {
                p = FoFNext[p];
                continue;
            }
            d = SphP[p].Density / Time3 / RhoBaryon;
            t = SphP[p].Temp;
            m = P[p].Mass;

            if ( t < 1e5 ) {
                if ( d>1e3 )
                    m_condensed += m;
                else
                    m_diffuse += m;
            }

            if ( t>=1e5 && t<1e7 )
                m_warmhot += m;

            if ( t>=1e7 )
                m_hot += m;

            p = FoFNext[p];
        }

        ep = group_pot( index );
#ifdef OUTPUTGROUP_DEBUG
        ek = group_ek( index );
#else
        ek = g->ek;
#endif

#ifdef OUTPUTGROUPLUM
        lum = group_luminosity( All.OutputGroupFreq, index, 1 );
#endif
        
        fprintf( fd, "%g,", g->mass );

        for( i=0; i<3; i++ )
            fprintf( fd, "%g,", g->cm[i] );

        fprintf( fd, "%g,%g,%g,%g,%g,%g,%g,",
        g->mass_table[0],
        g->mass_table[1],
        g->mass_table[4],
        m_diffuse,
        m_warmhot,
        m_hot,
        m_condensed
        );
#ifdef OUTPUTGROUPLUM
            fprintf( fd, "%g,", lum );
#endif

        fprintf( fd, "%g,%g,%g,%g,%g,%g\n",
                g->vr200, ek, ep, ek/ep, g->v_mean, g->v_disp );
#ifdef OUTPUTGROUP_DEBUG
        if ( index == 2 )
        break;
#endif

     }

    fclose( fd );
    put_end();
#ifdef OUTPUTGROUP_DEBUG
    endruns( "output_group-debug" );
#endif

}
#endif

void group_analysis() {

    int proj_tmp[3], i;
    char Sproj_tmp;

    put_header( "group analysis" );

#ifdef OUTPUTGROUP
    output_group();
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


void test_group_pot() {
#ifdef GROUP_POT_TEST
    double e, e_direct;
    int i;

    //tree_build();
    fof();
    for( i=0; i<10; i++ ) {

#ifdef GROUP_POT_TEST2
        i = 1;
#endif
        if ( ThisTask_Local == 0 )
            e = group_pot(i);

        MPI_Barrier( MpiComm_Local );
        e_direct = group_pot_direct(i);

        if ( ThisTask_Local == 0 )
        printf( "[%3i] pot fft: %g, pot direct: %g, err: %.2f%%, r: %g\n",
            i, e, e_direct, (e-e_direct)/e_direct*100, e/e_direct );
#ifdef GROUP_POT_TEST2
        break;
#endif
    }

    do_sync( "" );
    if ( ThisTask_Local == 0 )
        endruns( "group_pot-test" );

#endif

}
