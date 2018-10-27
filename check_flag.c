#include "allvars.h"

void check_flag() {

    if ( All.MF )
        if ( All.FoF == 0 )
            endrun( "FoF is required by MF Analysis!" )

    if ( All.Group ) {
        if ( All.FoF == 0 )
            endrun( "FoF is required by Group Analysis!" )
        check_group_flag();
    }


}
