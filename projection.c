#include "allvars.h"

void init_projection() {
    if ( !( All.ProjectDirection == 0 ||
            All.ProjectDirection == 1 ||
            All.ProjectDirection == 2 ) ) {
        printf( "project direction must be 0, 1 or 2 !\n" );
        endrun( 20171207 );
    }
    writelog( "initialize projection ...\n" );
    switch( All.ProjectDirection ) {
        case 0:
            proj_i = 1;
            proj_j = 2;
            proj_x = All.End[1] - All.Start[1];
            proj_y = All.End[2] - All.Start[2];
            break;
        case 1:
            proj_i = 0;
            proj_j = 2;
            proj_x = All.End[0] - All.Start[0];
            proj_y = All.End[2] - All.Start[2];
            break;
        case 2:
            proj_i = 0;
            proj_j = 1;
            proj_x = All.End[0] - All.Start[0];
            proj_y = All.End[1] - All.Start[1];
            break;
    }
    if ( proj_x != proj_y ) {
        printf( "proj_x != proj_y is no support !!!\n" );
        endrun( 20171207 );
    }
    proj_size = proj_x;
    writelog( "proj_i = %d, proj_j = %d\n"
                     "proj_x = %g, proj_y = %g\n"
                     "proj_size = %g\n",
                     proj_i, proj_j,
                     proj_x, proj_y,
                     proj_size );
    writelog( "initialize projection ... done.\n" );
    writelog( sep_str );
}

