#include "allvars.h"
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

#ifdef GROUPKIN
double group_kin( long index ) {

    long p;
    double e;
    struct group_properties *g;

    g = &Gprops[index];

    p = g->Head;
    e = 0;
    while( p >= 0 ) {
        e += 0.5 * P[p].Mass * 
        ( SQR( P[p].Vel[0] ) + 
                    SQR( P[p].Vel[1] ) + 
                    SQR( P[p].Vel[2] ) );
        p = FoFNext[p];

    }
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
     double v_disp, ek, ep;
     create_dir( "%sVirial", GroupDir  );
     fd = myfopen( "w", "%sVirial/Virial_%03i.dat",
        GroupDir, SnapIndex);
     for ( index=0; index<Ngroups; index++  ) {
         if ( !group_present( index  )  )
             break;
             v_disp = group_vel_disp( index );
             ep = group_pot( index );
             ek = group_kin( index );
             fprintf( fd, "%g %g %g\n", ep, ek, v_disp );
     }
     fclose( fd );
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


