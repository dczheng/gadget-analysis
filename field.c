#include"allvars.h"

void field_cren_T_dens() {

    double *ntd_x, *ntd_y, *ntd_z, mm[4];
    long p; 
    int flag;

    mymalloc1( ntd_x, sizeof(double) * N_Gas );
    mymalloc1( ntd_y, sizeof(double) * N_Gas );
    mymalloc1( ntd_z, sizeof(double) * N_Gas );
    for( p=0; p<N_Gas; p++ ) {
        ntd_x[p] = SphP[p].Density / All.RhoBaryon;
        ntd_y[p] = SphP[p].Temp;
        ntd_z[p] = SphP[p].CRE_n * SphP[p].Density / guc.m_e * ( 1/CUBE(g2c.cm) ) ;
    }


    flag = 0;
    flag |= PDF2D_BIT_XLOG;
    flag |= PDF2D_BIT_YLOG;

    if ( All.FieldCrenTDensDensMin > 0 || 
            All.FieldCrenTDensDensMax > 0 ) {
        if ( All.FieldCrenTDensDensMin * All.FieldCrenTDensDensMax == 0 )
            endrun( 20180810 );

        flag |= PDF2D_BIT_FIXEDX;
        mm[0] = All.FieldCrenTDensDensMin;
        mm[1] = All.FieldCrenTDensDensMax;

    }

    if ( All.FieldCrenTDensTMin > 0 || 
            All.FieldCrenTDensTMax > 0 ) {
        if ( All.FieldCrenTDensTMin * All.FieldCrenTDensTMax == 0 )
            endrun( 20180810 );

        flag |= PDF2D_BIT_FIXEDY;
        mm[2] = All.FieldCrenTDensTMin;
        mm[3] = All.FieldCrenTDensTMax;

    }

    field2d( ntd_x, ntd_y, ntd_z, N_Gas, "FiledCrenTDens", flag, mm, 1000  );

    myfree( ntd_x );
    myfree( ntd_y );
    myfree( ntd_z );
}

