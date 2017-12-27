#include "allvars.h"

long npart, last, parent, offset, father;
int pt;

void tree_allocate() {
    size_t bytes;
    MaxNodes = npart * All.TreeAllocFactor;
    sprintf( LogBuf, "particle `%i` npartber = %ld, TreeAllocateFactor = %g, "
            "MaxNodes = %ld", pt, npart, All.TreeAllocFactor, MaxNodes );
    print_log( LogBuf );
    print_log( "allocate memory for tree" );
    if ( !( Nodes_Base = malloc( bytes = (MaxNodes+1) * sizeof( struct NODE ) ) ) ) {
        printf( "failed to allocate memory for %ld tree-nodes `Nodes_Base`(%g MB)\n", MaxNodes, bytes / 1024.0 / 1024.0 );
        endrun( 20171225 );
    }
    sprintf( LogBuf, "allocate memory for %ld tree-nodes `Nodes_Base`( %g MB )", MaxNodes, bytes / 1024.0 / 1024.0 );
    Nodes = Nodes_Base - npart;
    print_log( LogBuf );

    if ( !( NextNode = malloc( bytes = (npart) * sizeof( long ) ) ) ) {
        printf( "falied to allocate memory for `NextNode` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
        endrun( 20171227 );
    }
    sprintf( LogBuf, "allocate memory for `NextNode` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );
    print_log( sep_str );
}

void tree_free() {
    print_log( "free memory for tree" );
    print_log( sep_str );
    free( Nodes_Base );
    free( NextNode );
}

void tree_build_single() {
    long i, j, subnode, bits, nfree, n, nn, ii;
    struct NODE *nfreep;
    double max[3], min[3], len, lenhalf;
    nfree = npart;
    nfreep = &Nodes[nfree];
    print_log( "tree build ..." );
    for ( j=0; j<3; j++ ) {
        max[j] = DBL_MIN;
        min[j] = DBL_MAX;
    }
    for ( i=offset; i<offset+npart; i++ ) {
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
    for ( i=offset+1; i<offset+npart; i++ ) {
        ii = i - offset;
        if ( !P[i].Flag ) continue;
        n = npart;
        while(1) {
            if ( n >= npart ){
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
                if ( nfree-npart >= MaxNodes ){
                    printf( "Max number of tree nodes reached.\n" );
                    endrun( 20171225 );
                }
            }
        }
    }
    print_log( "tree build ... done." );
}

void tree_walk_recursive( long n, long sib, long father ) {
    int i, j;
    long nextsib, p, pp;
    if ( n >= npart ) {
        if ( last >= 0 ) {
            if ( last >= npart )
                Nodes[last].nextnode = n;
            else
                NextNode[last] = n;
        }
        last = n;
        /*
        debug_l[0] = n;
        for ( i=0; i<8; i++ )
            printf( "%li ", Nodes[n].suns[i] );
        printf( "\n" );
        */
        for ( i=0; i<8; i++ ) {
            if ( (p = Nodes[n].suns[i]) >= 0 ) {
                for ( j=i+1; j<8; j++ )
                    if ( (pp = Nodes[n].suns[j]) >= 0 )
                        break;
                nextsib = ( j<8 ) ? pp : sib;
                tree_walk_recursive( p, nextsib, n );
                Nodes[n].sibling = sib;
                Nodes[n].father = father;
            }
        }
    }
    else {
        if ( last >= 0 ) {
            if ( last >= npart ) {
                Nodes[last].nextnode = n;
            }
            else {
                NextNode[last] = n;
            }
        }
        last = n;
    }
}

void tree_walk_test(){
    long n, signal;
    n = npart;
    signal = 0;
    while ( n>=0 ) {
        if ( n<npart ){
            printf( "%li -> ", n );
            if ( !(++signal % 10)  )
                printf( "\n" );
            n = NextNode[n];
        }
        else {
            n = Nodes[n].nextnode;
        }
    }
}

void tree_build( int ptype ) {
    pt = ptype;
    int i, j;
    long num;
    npart = get_particle_num( pt );
    if ( npart == 0 ) {
        printf( "particle `%li` npartber is zero !!!\n", pt );
        endrun( 20171225 );
    }
    /*
    npart = 30;
    All.TreeAllocFactor = 2;
    */
    offset = get_particle_offset( pt );
    for ( i=offset, num=0; i<offset+num; i++ )
        if ( P[i].Type ) num++;
    sprintf( LogBuf, "%i particle used in tree build", num );
    print_log( LogBuf );
    tree_allocate();
    tree_build_single();
    /*
    for ( i=0; i<npart*All.TreeAllocFactor; i++ ) {
        printf( "%2i %2i: ", i, i+npart );
        for ( j=0; j<8; j++ )
            printf( "%3li ", Nodes_Base[i].suns[j] );
        printf( "\n" );
    }
    */
    last = -1;
    print_log( "tree walk build ..." );
    tree_walk_recursive( npart, -1, -1 );
    if ( last >= npart ) {
        Nodes[last].nextnode = -1;
    }
    else {
        NextNode[last] = -1;
    }
    //tree_walk_test();
    print_log( "tree walk build ... done." );
    print_log( sep_str );
}
