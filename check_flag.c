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

    writelog( "check parameters ...\n" );
    A_NEED_B( All.FoF, All.Tree );
    A_NEED_B( All.MF, All.FoF);

    A_NEED_B( All.Group, All.FoF );

    A_NEED_B( All.BPdf, All.ReadB );
    A_NEED_B( All.GroupSpec, All.RadSpec );
    A_NEED_B( All.GroupU, All.Readu );

    A_NEED_B( All.RadSpec, All.ReadB );
    A_NEED_B( All.RadSpec, All.ReadCre );

    A_NEED_B( All.TotSpec, All.RadSpec );

    A_NEED_B( All.CrePressurePdf, All.ReadCre );
    A_NEED_B( All.CrePressurePdf, All.ReadCr );
    A_NEED_B( All.CrePressurePdf, All.Readu );

    A_NEED_B( All.RadSlice, All.RadSpec );

    if ( All.CorrTdiffDens && NTask != 2 ) {
        printf( "`NTask = 2` is required by `CorrTdiffDens`\n" );
        endrun( 20190720 );
    }

    if ( All.PdfTdiffDens && NTask != 2 ) {
        printf( "`NTask = 2` is required by `CorrTdiffDens`\n" );
        endrun( 20190730 );
    }

    if ( All.NGrid <= 1 ) {
        printf( "`NGrid > 1` ! is required\n" );
        endrun( 20190720 );
    }

    A_NEED_B( All.CorrTdiffDens, All.NGrid );
    A_NEED_B( All.PdfTdiffDens,  All.NGrid );
    A_NEED_B( All.CorrGas, All.NGrid );
    A_NEED_B( All.CorrDM, All.NGrid );

    A_NEED_B( All.TPdf, All.TPdfN );
    A_NEED_B( All.DensPdf, All.DensPdfN );
    A_NEED_B( All.RadSpec, All.ReadHsml );
    A_NEED_B( All.BSmooth, All.ReadHsml );
    A_NEED_B( All.BSmooth, All.ReadB );
    A_NEED_B( All.BSmooth, All.Tree );

    if ( All.TPdf && All.TPdfN < 1 ) {
        printf( "`TPdfN > 1` is required by TPdf\n" );
        endrun( 20190801 );
    }

    if ( All.DensPdf && All.DensPdfN < 1 ) {
        printf( "`DensPdfN > 1` is required by DensPdf\n" );
        endrun( 20190801 );
    }

    A_NEED_B( All.GroupTempStack, All.GroupTempStackRmin );
    A_NEED_B( All.GroupTempStack, All.GroupTempStackRmax );
    A_NEED_B( All.GroupTempStack, All.GroupTempStackRN );
    A_NEED_B( All.GroupTempStack, All.Group );
    A_NEED_B( All.GroupTempStack, All.FoF );
    A_NEED_B( All.GroupTempStack, All.Tree );

    A_NEED_B( All.GroupTempProfile, All.Group );
    A_NEED_B( All.GroupTempProfile, All.FoF );
    A_NEED_B( All.GroupTempProfile, All.Tree );
    A_NEED_B( All.GroupTempProfile, All.GroupTempProfileRmin );
    A_NEED_B( All.GroupTempProfile, All.GroupTempProfileRmax );

    A_NEED_B( All.CrenTPdf, All.ReadCre );
    A_NEED_B( All.BDensPdf, All.ReadB );
    A_NEED_B( All.FieldCrenTDens, All.ReadCre );
    A_NEED_B( All.DivBErrPdf, All.ReadB );
    A_NEED_B( All.DivBErrPdf, All.ReadDivB );
    A_NEED_B( All.DivBErrPdf, All.ReadHsml );


    A_NEED_B( All.GroupGasRatio, All.FoF );
    A_NEED_B( All.GroupGasRatio, All.Group );
    A_NEED_B( All.HsmlTPdf, All.ReadHsml );
    A_NEED_B( All.UTPdf, All.Readu );
    A_NEED_B( All.GroupCre, All.Readu );

    A_NEED_B( All.CREeSlice, All.Readu );
    A_NEED_B( All.CREeSlice, All.ReadCre );
    A_NEED_B( All.MachSlice, All.ReadMach );
    A_NEED_B( All.BSlice, All.ReadB );


}
