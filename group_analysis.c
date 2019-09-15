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

void group_analysis() {

    int proj_tmp[3], i;
    char Sproj_tmp;

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

}

#endif
