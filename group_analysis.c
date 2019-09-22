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
            fprintf( fd, "%e ", get_particle_radio( p, i ) );
        }
        fprintf( fd, "\n" );
        p = FoFNext[p];
    }
    fclose( fd );
    //hendruns( "output_group" );

}

void group_analysis() {

    int proj_tmp[3], i;
    char Sproj_tmp;

    put_header( "group analysis" );

    output_group( 0 );

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

#ifdef GROUPELECSPEC
        group_electron_spectrum();
        put_sep0;
#endif


#ifdef GROUPSPEC
        group_spectrum();
        //group_spectrum_index();
#endif

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
