#include "allvars.h"

void magnetic_field_analysis() {
     long i, index;
     double sum, b, bmax, min, max;
     sum = 0;
     index = 0;
     bmax = 0;
     for ( i=0; i<Particle[0].num; i++ ) {
         /*
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
                 */
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

void density_analysis() {
     long i, index;
     double sum, min, max, rho;
     fputs( sep_str, stdout );
     fputs( "analyze density ...\n", stdout );
     sum = 0;
     index = 0;
     /*
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
     min = 1e10;
     max = -1;
     for ( i=0; i<Particle[0].num; i++ ) {
        min = ( Particle[0].rho[i]<min ) ? Particle[0].rho[i] : min;
        max = ( Particle[0].rho[i]>max ) ? Particle[0].rho[i] : max;
     }
     fprintf( stdout, "min = %e, max = %e\n", min, max );
     */
    sum = 0;
     for ( i=0; i<Particle[0].num; i++ ) {
         sum += Particle[0].m[i];
     }
     fprintf( stdout, "m=%e\n", sum );
     fputs( sep_str, stdout );
}

int compare_for_sort_group_by_mass( const void *a, const void *b ) {
    return (( ( struct group_struct* )a )->Mass < ( ( struct group_struct* )b )->Mass ) ? 1 : -1;
}

void sort_group_by_mass() {
    qsort( (void*)group, TotNgroups, sizeof( struct group_struct ), &compare_for_sort_group_by_mass );
}

void show_group_info() {
    int i;
    for ( i=0; i<TotNgroups; i++ ) {
        if ( group[i].Mass >= 1e3 )
        fprintf( stdout, "mass: %10.3f, cm: ( %10.3f, %10.3f, %10.3f ), vel: ( %10.3f, %10.3f, %10.3f )\n",
                group[i].Mass, group[i].CM[0], group[i].CM[1], group[i].CM[2],
                group[i].Vel[0], group[i].Vel[1], group[i].Vel[2] );
    }
}

void group_analysis() {
    read_group();
    sort_group_by_mass();
    show_group_info();
    free_group();
}

void velocity_analysis() {
    int pt;
    long i;
    float vmax[3], v;
    pt = 0;
    vmax[0] = Particle[pt].vel[0];
    vmax[1] = Particle[pt].vel[1];
    vmax[2] = Particle[pt].vel[2];
    v = sqrt( pow( vmax[0], 2 ) + pow( vmax[1], 2 ) + pow( vmax[2], 2 ) );
    for ( i=0; i<Particle[pt].num; i++ ) {
        if ( sqrt( pow( Particle[pt].vel[ i*3+0 ], 2 ) +
                   pow( Particle[pt].vel[ i*3+1 ], 2 ) +
                   pow( Particle[pt].vel[ i*3+2 ], 2 ) ) > v ) {
            vmax[0] = Particle[pt].vel[ i*3+0 ];
            vmax[1] = Particle[pt].vel[ i*3+1 ];
            vmax[2] = Particle[pt].vel[ i*3+2 ];
            v = sqrt( pow( vmax[0], 2 ) +
                    pow( vmax[1], 2 ) +
                    pow( vmax[2], 2 ) );
        }
    }
    fprintf( stdout, "Vmax: %f ( %f, %f, %f )\n", v, vmax[0],
            vmax[1], vmax[2] );
}

