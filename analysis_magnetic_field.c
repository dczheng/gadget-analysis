#include "allvars.h"

void Analysis_Magnetic_Field() {
     long i;
     double sum;
     sum = 0;
     for ( i=0; i<Particle[0].num; i++ ) {
         sum += pow( Particle[0].mag[i*3+0], 2 ) +
                pow( Particle[0].mag[i*3+1], 2 ) +
                pow( Particle[0].mag[i*3+2], 2 );
     }

     fprintf( stdout, "Total Magnetic Field: %lf\n", sqrt( sum / Particle[0].num ) );
}
