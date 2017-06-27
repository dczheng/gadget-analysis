#include "allvars.h"

void magnetic_field_analysis() {
     long i, index;
     double sum, b, bmax, min, max;
     sum = 0;
     index = 0;
     bmax = 0;
     for ( i=0; i<Particle[0].num; i++ ) {
         if( Particle[0].mag[i*3+0] != 0 ||
             Particle[0].mag[i*3+1] != 0 ||
             Particle[0].mag[i*3+2] != 0 )
         fprintf( stdout, "( %.2f, %.2f, %.2f ): %e %e %e\n",
                 Particle[0].pos[i*3+0],
                 Particle[0].pos[i*3+1],
                 Particle[0].pos[i*3+2],
                 Particle[0].mag[i*3+0],
                 Particle[0].mag[i*3+1],
                 Particle[0].mag[i*3+2] );
         b = pow( Particle[0].mag[i*3+0], 2 ) +
             pow( Particle[0].mag[i*3+1], 2 ) +
             pow( Particle[0].mag[i*3+2], 2 );
         if ( b>bmax ) {
             index =i;
             bmax = b;
         }
         sum += b;
     }
    fprintf( stdout, "max magnetic field: ( %.2f, %.2f, %.2f ): %e %e %e %e\n",
            Particle[0].pos[index*3+0],
            Particle[0].pos[index*3+1],
            Particle[0].pos[index*3+2],
            Particle[0].mag[index*3+0],
            Particle[0].mag[index*3+1],
            Particle[0].mag[index*3+2], sqrt( bmax ) );
     fprintf( stdout, "Total Magnetic Field: %e\n", sqrt( sum / Particle[0].num ) );
}
