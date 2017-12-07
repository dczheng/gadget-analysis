#include "allvars.h"

void slice() {
    char buf[200];
    int pt;
    long long offset, num, index, i;
    struct particle_data p_tmp;
    struct sph_particle_data sphp_tmp;
    print_log( "determine slice info ..." );

    if ( !( para.ProjectDirection == 1 ||
            para.ProjectDirection == 2 ||
            para.ProjectDirection == 3 ) ) {
        printf( "project direction must be 1, 2 or 3 !\n" );
        endrun( 20171207 );
    }

    if ( para.EndX == 0 ){
        para.StartX = 0;
        para.EndX = BoxSize;
    }

    if ( para.EndY == 0 ){
        para.StartY = 0;
        para.EndY = BoxSize;
    }

    if ( para.EndZ == 0 ){
        para.StartZ = 0;
        para.EndZ = BoxSize;
    }
    sprintf( LogBuf, "StartX: %g, EndX: %g\n"
            "StartY: %g, EndY: %g\n"
            "StartZ: %g, EndZ: %g",
            para.StartX, para.EndX,
            para.StartY, para.EndY,
            para.StartZ, para.EndZ );
    print_log( LogBuf );

    for ( pt=0,offset=0; pt<6; pt++ ) {
        index = offset;
        num = header.npartTotal[pt] + ( ( (long long)header.npartTotalHighWord[pt] ) << 32 );
        if ( num == 0 ) {
            SliceStart[pt] = -1;
            SliceEnd[pt] = -1;
            continue;
        }
        for ( i=offset; i<offset+num; i++ ) {
            if ( P[i].Pos[0] > para.StartX &&
                 P[i].Pos[0] < para.EndX &&
                 P[i].Pos[1] > para.StartY &&
                 P[i].Pos[1] < para.EndY &&
                 P[i].Pos[2] > para.StartZ &&
                 P[i].Pos[2] < para.EndZ ) {
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
        SliceEnd[pt] = offset + index - 1;
        offset += num;
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

void get_projection_index( int *i, int *j ) {
    switch( para.ProjectDirection ) {
        case 1:
            *i = 2;
            *j = 3;
            break;
        case 2:
            *i = 1;
            *j = 3;
            break;
        case 3:
            *i = 1;
            *j = 2;
            break;
    }
}

