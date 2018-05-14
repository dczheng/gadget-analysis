#include "allvars.h"

#define FACT1 0.366025403785	/* FACT1 = 0.5 * (sqrt(3)-1) */
#define FACT2 0.86602540        /* FACT2 = 0.5 * sqrt(3) */

int ngb_fof( double *searchcenter, double h ) {
    long n, p;
    int ngbnum;
    double dx, dy, dz, dist, r2;
    struct NODE *c;
    ngbnum = 0;
    n = NumPart;
    while ( n>=0 ) {
        if ( n < NumPart ) {
            if ( !( ( 1 << P[n].Type ) & All.TreePartType ) )
                continue;

            p = n;
            n = NextNode[n];
            dist = h;

            dx = NGB_PERIODIC( P[p].Pos[0] - searchcenter[0] );
            if ( dx > dist ) continue;

            dy = NGB_PERIODIC( P[p].Pos[1] - searchcenter[1] );
            if ( dy > dist ) continue;

            dz = NGB_PERIODIC( P[p].Pos[2] - searchcenter[2] );
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

            dist += FACT1 * c->len;
            r2 = dx*dx + dy*dy + dz*dz;
            if ( r2 > dist*dist )
                continue;

            dist = h - FACT2 * c -> len;
            if ( dist > 0 )
                if ( r2 < dist * dist ) {
                    if ( c -> bitflags & ( 1 << BITFLAG_INSIDE_LINKINGLENGTH ) ) {
                        p = c -> nextnode;
                        while( p >=0 ) {
                            if ( p < NumPart ) {
                                if ( ( 1 << P[p].Type ) & ( All.TreePartType ) ) {
                                    dx = NGB_PERIODIC( P[p].Pos[0] - searchcenter[0] );
                                    dy = NGB_PERIODIC( P[p].Pos[1] - searchcenter[1] );
                                    dz = NGB_PERIODIC( P[p].Pos[2] - searchcenter[2] );

                                    if ( dx*dx + dy*dy + dz*dz > h*h )
                                        break;

                                    Ngblist[ngbnum++] = p;
                                    break;
                                }
                                p = NextNode[p];
                            }
                            else
                                p = Nodes[p].nextnode;
                        }
                        continue;
                    }
                    else {
                        c -> bitflags |= ( 1 << BITFLAG_INSIDE_LINKINGLENGTH );
                    }
                }
            n = c->nextnode;
        }
    }
    return ngbnum;
}

int ngb( double *searchcenter, double h ) {
    long n, p;
    int ngbnum;
    double dx, dy, dz, dist;
    struct NODE *c;
    ngbnum = 0;
    n = NumPart;
    while ( n>=0 ) {
        if ( n < NumPart ) {
            if ( !( ( 1 << P[n].Type ) & All.TreePartType ) )
                continue;
            p = n;
            n = NextNode[n];
            dist = h;

            dx = NGB_PERIODIC( P[p].Pos[0] - searchcenter[0] );
            if ( dx > dist ) continue;

            dy = NGB_PERIODIC( P[p].Pos[1] - searchcenter[1] );
            if ( dy > dist ) continue;

            dz = NGB_PERIODIC( P[p].Pos[2] - searchcenter[2] );
            if ( dz > dist ) continue;

            if ( dx*dx + dy*dy + dz*dz >  dist*dist )
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

            dist += FACT1 * c->len;
            if ( dx*dx + dy*dy + dz*dz > dist*dist )
                continue;

            n = c->nextnode;
        }
    }
    return ngbnum;
}
