#include "allvars.h"

long npart, last, parent, offset, father;
int pt;

void tree_allocate() {
    size_t bytes;
    MaxNodes = npart * All.TreeAllocFactor;
    writelog( "particle `%i` npartber = %ld, TreeAllocateFactor = %g, "
            "MaxNodes = %ld\n", pt, npart, All.TreeAllocFactor, MaxNodes );
    writelog( "allocate memory for tree\n" );
    mymalloc( Nodes_Base, ( MaxNodes+1 ) * sizeof( struct NODE ) );
    Nodes = Nodes_Base - npart;

    mymalloc( NextNode, (npart) * sizeof( long ) );
    writelog( sep_str );
}

void tree_free() {
    writelog( "free memory for tree\n" );
    myfree( Nodes_Base );
    myfree( NextNode );
    writelog( sep_str );
}

void tree_build_single() {
    long i, j, subnode, bits, nfree, n, nn, ii;
    struct NODE *nfreep;
    double max[3], min[3], len, lenhalf;
    nfree = npart;
    nfreep = &Nodes[nfree];
    writelog( "tree build ...\n" );
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
    for ( i=0; i<3; i++ ) {
        writelog( "max[%i]=%g, min[%i]=%g\n",
                i, max[i], i, min[i] );
    }
    len *= 1.001;
    writelog( "len=%g\n", len );
    /* initialize first node */
    for ( j=0; j<3; j++ ) {
        nfreep->center[j] = ( max[j] + min[j] ) * 0.5;
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
                    subnode += (( P[n+offset].Pos[j] > nfreep->center[j] ) ? bits : 0);
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
    writelog( "tree build ... done.\n" );
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
    FILE *fd;
    long n, signal;
    n = npart;
    signal = 0;
    fd = fopen( "walk.txt", "w" );
    while ( n>=0 ) {
        if ( n<npart ){
            fprintf( fd, "%li\n", n );
            n = NextNode[n];
        }
        else {
            n = Nodes[n].nextnode;
        }
    }
    fclose( fd );
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
    writelog( "%i particle used in tree build\n", num );
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
    writelog( "tree walk build ...\n" );
    tree_walk_recursive( npart, -1, -1 );
    if ( last >= npart ) {
        Nodes[last].nextnode = -1;
    }
    else {
        NextNode[last] = -1;
    }
    //tree_walk_test();
    writelog( "tree walk build ... done.\n" );
    writelog( sep_str );
}
