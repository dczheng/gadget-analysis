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
    A_NEED_B( All.GroupSpec, All.RadSpec );

    if ( All.Group )
        check_group_flag();

    A_NEED_B( All.RadSpec, All.ReadB );
    A_NEED_B( All.RadSpec, All.ReadCre );

    A_NEED_B( All.TotSpec, All.RadSpec );

    A_NEED_B( All.CrePressurePdf, All.ReadCre );
    A_NEED_B( All.CrePressurePdf, All.ReadCr );
    A_NEED_B( All.CrePressurePdf, All.Readu );

    A_NEED_B( All.Phase, All.ReadElec );
    A_NEED_B( All.Phase, All.Readu );

    A_NEED_B( All.RadSlice, All.RadSpec );

    A_NEED_B( All.GroupTemp, All.GroupTempBins );
    A_NEED_B( All.GroupTemp, All.GroupTempRmin );

    if ( All.CorrTdiffDens && NTask != 2 ) {
        printf( "`Master = 2` is required by All.CorrTdiffDens\n" );
        endrun( 20190720 );
    }

    if ( All.CorrGrid <= 1 ) {
        printf( "All.CorrGrid <= 1 !\n" );
        endrun( 20190720 );
    }

    A_NEED_B( All.CorrTdiffDens, All.CorrGrid );

    A_NEED_B( All.CorrGrid, All.CorrRN );

}
