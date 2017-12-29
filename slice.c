#include "allvars.h"

void slice() {
    char buf[200];
    int pt;
    long long offset, num, index, i;
    struct particle_data p_tmp;
    struct sph_particle_data sphp_tmp;
    print_log( "determine slice info ..." );

    if ( All.End[0] == 0 ){
        All.Start[0] = 0;
        All.End[0] = All.BoxSize;
    }

    if ( All.End[1] == 0 ){
        All.Start[1] = 0;
        All.End[1] = All.BoxSize;
    }

    if ( All.End[2] == 0 ){
        All.Start[2] = 0;
        All.End[2] = All.BoxSize;
    }
    sprintf( LogBuf, "StartX: %g, EndX: %g\n"
            "StartY: %g, EndY: %g\n"
            "StartZ: %g, EndZ: %g",
            All.Start[0], All.End[0],
            All.Start[1], All.End[1],
            All.Start[2], All.End[2] );
    print_log( LogBuf );

    for ( pt=0; pt<6; pt++ ) {
        offset = get_particle_offset( pt );
        num = get_particle_num( pt );
        index = offset;
        if ( num == 0 ) {
            SliceStart[pt] = -1;
            SliceEnd[pt] = -1;
            continue;
        }
        sprintf( LogBuf, "particle %i: offset=%li, num=%li",
                pt, offset, num );
        print_log( LogBuf );
        for ( i=offset; i<offset+num; i++ ) {
            if ( P[i].Pos[0] >= All.Start[0] &&
                 P[i].Pos[0] <= All.End[0] &&
                 P[i].Pos[1] >= All.Start[1] &&
                 P[i].Pos[1] <= All.End[1] &&
                 P[i].Pos[2] >= All.Start[2] &&
                 P[i].Pos[2] <= All.End[2] ) {
                p_tmp = P[index];
                P[index] = P[i];
                P[i] = P[index];
                if ( pt == 0 ) {
                    sphp_tmp = SphP[ index ];
                    SphP[index] = SphP[i];
                    SphP[i] = sphp_tmp;
                }
                index ++;
            }
        }
        SliceStart[pt] = offset;
        SliceEnd[pt] = index;
    }
    sprintf( LogBuf, "Slice Start: " );
    for ( pt=0; pt<6; pt++ ) {
        sprintf( LogBuf, "%s%ld ", LogBuf, SliceStart[pt] );
    }
    print_log( LogBuf );

    sprintf( LogBuf, "Slice End: " );
    for ( pt=0; pt<6; pt++ ) {
        sprintf( LogBuf, "%s%ld ", LogBuf, SliceEnd[pt] );
    }
    print_log( LogBuf );

    print_log( "determine slice info ... done." );
    print_log( sep_str );
}

