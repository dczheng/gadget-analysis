#include "allvars.h"

#define A_NEED_B( A, B ){\
    if ( A ) { \
        if ( B == 0 ) { \
            printf( "`%s` is required by `%s` !\n", #B, #A );\
            endrun( 20181031 );\
        } \
    } \
}

void check_flag() {

    A_NEED_B( All.MF, All.FoF);

    A_NEED_B( All.Group, All.FoF );

    A_NEED_B( All.BPdf, All.ReadB );

    if ( All.Group )
        check_group_flag();

    A_NEED_B( All.RadSpec, All.ReadB );
    A_NEED_B( All.RadSpec, All.ReadHge );

    A_NEED_B( All.TotSpec, All.RadSpec );

    A_NEED_B( All.HgePressurePdf, All.ReadHge );
    A_NEED_B( All.HgePressurePdf, All.ReadCr );


}
