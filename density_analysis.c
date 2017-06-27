#include "allvars.h"

void density_analysis() {
     long i, index;
     double sum, min, max, rho;
     fputs( sep_str, stdout );
     fputs( "analyze density ...\n", stdout );
     sum = 0;
     index = 0;
     for ( i=0; i<Particle[0].num; i++ ) {
         rho = Particle[0].rho[i] * 1.989 * 1e43 / pow( 3.08e21,3 ); //* 150000 * 3.085e21;
         fprintf( stdout, "%e\n", rho );
         if ( i==10 ) break;
     }
     max = -1;
     min = 1000;
     for ( i=0; i<Particle[0].num; i++ ) {
        min = ( Particle[0].m[i]<min ) ? Particle[0].m[i] : min;
        max = ( Particle[0].m[i]>max ) ? Particle[0].m[i] : max;
     }
     fprintf( stdout, "min = %lf, max = %lf\n", min, max );
     fputs( sep_str, stdout );
}
