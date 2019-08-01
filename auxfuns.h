
#define get_B( i ) ( sqrt(SQR(SphP[i].B[0]) + SQR(SphP[i].B[1]) + SQR(SphP[i].B[2]) ) )
#define get_pressure( i ) (SphP[i].u * ( GAMMA_MINUS1 * SphP[i].Density ))

#define get_gas_density_min_max( dmin, dmax ) {\
    long p;\
    dmin = DBL_MAX;\
    dmax = -dmin;\
    for( p=0; p<N_Gas; p++ ) {\
        vmax2( dmax, SphP[p].Density );\
        vmin2( dmin, SphP[p].Density );\
        }\
}

#define get_gas_temp_min_max( Tmin, Tmax ) {\
    long p;\
    Tmin = DBL_MAX;\
    Tmax = -Tmin;\
    for( p=0; p<N_Gas; p++ ) {\
        vmax2( Tmax, SphP[p].Temp );\
        vmin2( Tmin, SphP[p].Temp );\
        }\
}
