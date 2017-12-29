#include "allvars.h"

double fabsx;
#define NGB_PERIODIC( x ) ( fabsx=fabs(x), ( (fabsx > All.HalfBoxSize) ? ( All.BoxSize-fabsx ) : fabsx ) )
#define FAC  0.366025403785  /* 0.5 * ( sqrt(3) - 1 ) */

long offset, npart;

int ngb_fof( double *searchcenter, double h, int pt ) {
    long n, p;
    int ngbnum;
    double dx, dy, dz, dist, h2;
    struct NODE *c;
    npart = get_particle_num( pt );
    offset = get_particle_offset( pt );
    ngbnum = 0;
    n = npart;
    h2 = h*h;
    while ( n>=0 ) {
        if ( n < npart ) {
            p = n;
            n = NextNode[n];
            dist = h;
            dx = NGB_PERIODIC( P[p+offset].Pos[0] - searchcenter[0] );
            if ( dx > dist ) continue;
            dy = NGB_PERIODIC( P[p+offset].Pos[1] - searchcenter[1] );
            if ( dy > dist ) continue;
            dz = NGB_PERIODIC( P[p+offset].Pos[2] - searchcenter[2] );
            if ( dz > dist ) continue;
            if ( dx*dx + dy*dy + dz*dz > dist*dist )
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
            if ( dx*dx + dy*dy + dz*dz > dist*dist )
                continue;
            n = c->nextnode;
        }
    }
    return ngbnum;
}
