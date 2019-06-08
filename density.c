#include "allvars.h"

void density_slice() {
    int num, i;
    char buf[100];

    writelog( "gas density silce ...\n" );
    num = All.SliceEnd[0] - All.SliceStart[0];
    mymalloc2( image.data, sizeof( double ) * num );

    for ( i=All.SliceStart[0]; i<num; i++ ) {
        image.data[i] = SphP[i].Density * ( g2c.g / CUBE( g2c.cm ) );
    }

    sprintf( buf, "%sDensity", All.OutputPrefix );
    create_dir( buf );
    sprintf( buf, "%s/Density_%.2f.dat", buf, All.RedShift );

    make_slice_img( 0 );

    write_img2( buf, "gas density slice" );
    myfree( image.data );

    writelog( "gas density silce ... done.\n" );
    put_sep;
}

int compare_gas_rho( const void *a, const void *b ){
    return ((((struct Sph_Particle_Data* )a)->Density) < (((struct Sph_Particle_Data*)b)->Density)) ? 1: -1;
}

void sort_gas_rho(){
    int i;
    qsort( (void*)SphP, N_Gas, sizeof( struct Sph_Particle_Data ), &compare_gas_rho );

    for ( i=0; i<10; i++ ) {
        printf( "%g\n", SphP[i].Density * ( g2c.g / CUBE(g2c.cm) ) );
    }

    printf( "\n" );

    for ( i=0; i<10; i++ ) {
        printf( "%g\n", SphP[N_Gas-10+i].Density * ( g2c.g / CUBE(g2c.cm) ) );
    }
}

