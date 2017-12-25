#include "allvars.h"

void init_projection() {
    if ( !( para.ProjectDirection == 0 ||
            para.ProjectDirection == 1 ||
            para.ProjectDirection == 2 ) ) {
        printf( "project direction must be 0, 1 or 2 !\n" );
        endrun( 20171207 );
    }
    print_log( "initialize projection ..." );
    switch( para.ProjectDirection ) {
        case 0:
            proj_i = 1;
            proj_j = 2;
            proj_x = para.End[1] - para.Start[1];
            proj_y = para.End[2] - para.Start[2];
            break;
        case 1:
            proj_i = 0;
            proj_j = 2;
            proj_x = para.End[0] - para.Start[0];
            proj_y = para.End[2] - para.Start[2];
            break;
        case 2:
            proj_i = 0;
            proj_j = 1;
            proj_x = para.End[0] - para.Start[0];
            proj_y = para.End[1] - para.Start[1];
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

