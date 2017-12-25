#include "allvars.h"

void tree_allocate( int pt ) {
    long num;
    size_t bytes;
    num = get_particle_num( pt );
    MaxNodes = num * All.TreeAllocFactor;
    sprintf( LogBuf, "particle `%i` number = %ld, TreeAllocateFactor = %g, "
            "MaxNodes = %ld", pt, num, All.TreeAllocFactor, MaxNodes );
    print_log( LogBuf );
    print_log( "allocate memory for tree" );
    if ( !( Nodes_Base = malloc( bytes = (MaxNodes+1) * sizeof( struct NODE ) ) ) ) {
        printf( "failed to allocate memory for %ld tree-nodes (%g MB)\n", MaxNodes, bytes / 1024.0 / 1024.0 );
        endrun( 20171225 );
    }
    sprintf( LogBuf, "allocate memory for %ld tree-nodes ( %g MB )", MaxNodes, bytes / 1024.0 / 1024.0 );
    Nodes = Nodes_Base - num;
    print_log( LogBuf );
    print_log( sep_str );
}

void tree_free() {
    print_log( "free memory for tree" );
    print_log( sep_str );
    free( Nodes_Base );
}

void tree_build( int pt ) {
    long offset, num, i, j, subnode, bits, nfree, n, nn, parent, ii;
    struct NODE *nfreep;
    double max[3], min[3], len, lenhalf;
    tree_allocate( pt );
    offset = get_particle_offset( pt );
    num = get_particle_num( pt );
    if ( num == 0 ) {
        printf( "particle `%li` number is zero !!!\n", pt );
        endrun( 20171225 );
    }
    nfree = num;
    nfreep = &Nodes[nfree];
    print_log( "tree build ..." );
    for ( j=0; j<3; j++ ) {
        max[j] = DBL_MIN;
        min[j] = DBL_MAX;
    }
    for ( i=offset; i<offset+num; i++ ) {
        for ( j=0; j<3; j++ ) {
            max[j] = ( P[i].Pos[j] > max[j] ) ? P[i].Pos[j] : max[j];
            min[j] = ( P[i].Pos[j] < min[j] ) ? P[i].Pos[j] : min[j];
        }
    }
    len = DBL_MIN;
    for ( i=0; i<3; i++ ) {
        len = ( max[i] - min[i] > len ) ? ( max[i] - min[i] ) : len;
    }
    sprintf( LogBuf, "" );
    for ( i=0; i<3; i++ ) {
        sprintf( LogBuf, "%smax[%i]=%g, min[%i]=%g\n",
                LogBuf, i, max[i], i, min[i] );
    }
    len *= 1.01;
    sprintf( LogBuf, "%slen=%g", LogBuf, len );
    print_log( LogBuf );
    /* initialize first node */
    for ( j=0; j<3; j++ ) {
        nfreep->center[j] = ( max[j] + min[j] ) / 2;
    }
    for ( j=0; j<8; j++ ) {
        nfreep->suns[j] = -1;
    }
    nfreep->len = len;
    for ( j=0, subnode=0, bits=1; j<3; j++, bits<<=1 )
        subnode += ( P[offset].Pos[j] > nfreep->center[j] ) ? bits : 0;
    nfreep->suns[subnode] = 0;
    /* initialize first node*/
    nfree++;
    nfreep++;
    for ( i=offset+1; i<offset+num; i++ ) {
        ii = i - offset;
        n = num;
        while(1) {
            if ( n >= num ){
                for( j=0, subnode=0, bits=1; j<3; j++, bits<<=1 )
                    subnode += ( P[i].Pos[j] > Nodes[n].center[j] ) ? bits : 0;
                nn = Nodes[n].suns[subnode];
                if ( nn>=0 ) {
                    parent = n;
                    n = nn;
                }
                else {
                    Nodes[n].suns[subnode] = ii;
                    break;
                }
            }
            else {
                Nodes[parent].suns[subnode] = nfree;
                nfreep->len = 0.5 * Nodes[parent].len;
                lenhalf = 0.25 * Nodes[parent].len;
                for ( j=0, bits=1; j<3; j++, bits<<=1 )
                    nfreep->center[j] = Nodes[parent].center[j] + ( (subnode & bits) ? ( lenhalf ) : ( -lenhalf ) );
                for ( j=0; j<8; j++ )
                    nfreep->suns[j] = -1;
                for ( j=0, subnode=0, bits=1; j<3; j++, bits<<=1 )
                    subnode += (( P[n].Pos[j] > nfreep->center[j] ) ? bits : 0);
                nfreep->suns[subnode] = n;
                n = nfree;
                nfree++;
                nfreep++;
                if ( nfree-num >= MaxNodes ){
                    printf( "Max number of tree nodes reached.\n" );
                    endrun( 20171225 );
                }
            }
        }
    }

    print_log( "tree build ... done." );
    print_log( sep_str );
}
