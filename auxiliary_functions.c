#include "allvars.h"

long get_particle_num( int pt ) {
    return ( header.npartTotal[pt] + ( ( (long) header.npartTotalHighWord[pt] ) << 32 ) );
}

long get_particle_offset( int pt ) {
    long offset, i;
    for ( i=0, offset=0; i<pt; i++ ){
        offset += header.npartTotal[i];
        offset += ( ( (long)header.npartTotalHighWord[i] ) << 32 );
    }
    return offset;
}

double get_B( long i ) {
    return sqrt( SQR(SphP[i].B[0]) + SQR(SphP[i].B[1]) + SQR(SphP[i].B[2]) );
}

void init_sep_str() {
    memset( sep_str, '-', SEP_LEN-2 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';
}

