#include "allvars.h"

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

    if ( All.GroupEleSpec )
        group_electron_spectrum();

    put_sep0;

    if ( All.GroupSpec ) {
        group_spectrum();
        //group_spectrum_index();
    }

    if ( All.GroupTempProfile )
        group_temp_profile();

    if ( All.GroupTempStack )
        group_temp_stack();

    if ( All.GroupGasRatio )
        group_gas_ratio();

    reset_img();
    put_sep0;

}

