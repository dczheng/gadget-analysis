#include "allvars.h"

double fabsx;
#define NGB_PERIODIC( x ) ( fabsx=fabs(x), ( (fabsx > HalfBoxSize) ? ( BoxSize-fabsx ) : fabsx ) )
#define FAC  0.366025403785  /* 0.5 * ( sqrt(3) - 1 ) */

int ngb_fof( double searchcenter[3], double h, long target, long npart ) {
    long n, p;
    int ngbnum;
    double dx, dy, dz, dist, h2;
    struct NODE *c;
    ngbnum = 0;
    n = npart;
    h2 = h*h;
    while ( n>=0 ) {
        if ( n < npart ) {
            p = n;
            n = NextNode[n];
            dist = h;
            dx = NGB_PERIODIC( P[p].Pos[0] - searchcenter[0] );
            if ( dx > dist ) continue;
            dy = NGB_PERIODIC( P[p].Pos[1] - searchcenter[1] );
            if ( dy > dist ) continue;
            dz = NGB_PERIODIC( P[p].Pos[2] - searchcenter[2] );
            if ( dz > dist ) continue;
            if ( dx*dx + dy*dy + dz*dz > h2 )
                continue;
            Ngblist[ngbnum++] = p;
        }
        else {
            c = &Nodes[n];
            n = c->sibling;
            dist = h + 0.5 * c->len;
            dx = NGB_PERIODIC( c->center[0] - searchcenter[0] );
            if ( dx > dist ) continue;
            dy = NGB_PERIODIC( c->center[1] - searchcenter[1] );
            if ( dy > dist ) continue;
            dz = NGB_PERIODIC( c->center[2] - searchcenter[2] );
            if ( dz > dist ) continue;
            dist += FAC * c->len;
            if ( dx*dx + dy*dy + dz*dz > h2 )
                continue;
            n = c->nextnode;
        }
    }
    return ngbnum;
}
