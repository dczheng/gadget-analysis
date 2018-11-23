
#define get_B( i ) ( sqrt(SQR(SphP[i].B[0]) + SQR(SphP[i].B[1]) + SQR(SphP[i].B[2]) ) )

#define get_pressure( i ) (SphP[i].u * ( GAMMA_MINUS1 * SphP[i].Density ))
