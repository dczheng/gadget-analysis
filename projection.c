#include "allvars.h"

void init_projection() {
    if ( !( para.ProjectDirection == 1 ||
            para.ProjectDirection == 2 ||
            para.ProjectDirection == 3 ) ) {
        printf( "project direction must be 1, 2 or 3 !\n" );
        endrun( 20171207 );
    }
    print_log( "initialize projection ..." );
    switch( para.ProjectDirection ) {
        case 1:
            proj_i = 1;
            proj_j = 2;
            proj_x = para.EndY - para.StartY;
            proj_y = para.EndZ - para.StartZ;
            break;
        case 2:
            proj_i = 0;
            proj_j = 2;
            proj_x = para.EndX - para.StartX;
            proj_y = para.EndZ - para.StartZ;
            break;
        case 3:
            proj_i = 0;
            proj_j = 1;
            proj_x = para.EndX - para.StartX;
            proj_y = para.EndY - para.StartY;
            break;
    }
    if ( proj_x != proj_y ) {
        printf( "proj_x != proj_y is no support !!!\n" );
        endrun( 20171207 );
    }
    proj_size = proj_x;
    sprintf( LogBuf, "proj_i = %d, proj_j = %d, "
                     "proj_x = %g, proj_y = %g\n"
                     "proj_size = %g",
                     proj_i, proj_j,
                     proj_x, proj_y,
                     proj_size );
    print_log( LogBuf );
    print_log( "initialize projection ... done." );
    print_log( sep_str );
}

