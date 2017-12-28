#include "allvars.h"

long npart, offset;
double LinkL, rhodm;

void fof_allocate() {
    size_t bytes;
    if ( !(Head = malloc( bytes = npart * sizeof( long ) ) ) )
        printf( "failed to allocate memory for `Head` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `Head` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(Tail = malloc( bytes = npart * sizeof( long ) ) ) )
        printf( "failed to allocate memory for `Tail` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `Tail` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(Len = malloc( bytes = npart * sizeof( long ) ) ) )
        printf( "failed to allocate memory for `Len` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `Len` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(Next = malloc( bytes = npart * sizeof( long ) ) ) )
        printf( "failed to allocate memory for `Next` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `Next` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(Ngblist = malloc( bytes = npart * sizeof( long ) ) ) )
        printf( "failed to allocate memory for `Ngblist` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `Ngblist` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    print_log( sep_str );
}

void fof_free() {
    free( Head );
    free( Tail );
    free( Len );
    free( Next );
    free( Ngblist );
}

void fof_find_groups() {
    int ngbnum, k;
    double pos[3];
    long i, ii, p, s, j, ss;
    print_log( "start fof find groups ..." );
    for ( i=offset; i<offset+npart; i++ ) {
        ii = i-offset;
        if ( !P[i].Flag ) continue;
        for( j=0; j<3; j++ )
            pos[j] = P[i].Pos[j];
        ngbnum = ngb_fof( pos, LinkL, npart );
        for ( k=0; k<ngbnum; k++ ) {
            j = Ngblist[k];
            if ( Head[ii] != Head[j] ) {
                if ( Len[Head[ii]] > Len[Head[j]] ) {
                    p = ii;
                    s = j;
                }
                else{
                    p = j;
                    s = ii;
                }
                Next[Tail[Head[p]]] = Head[s];
                Tail[Head[p]] = Tail[Head[s]];
                Len[Head[p]] += Len[Head[s]];
                ss = Head[s];
                do
                    Head[ss] = Head[p];
                while( (ss=Next[ss]) >=0 );
            }
        }
    }
    print_log( "fof find groups ... done" );
}

void fof( int pt ) {
    long i;
    double masstot, mass;
    print_log( "start fof ..." );
    npart = get_particle_num( pt );
    offset = get_particle_offset( pt );
    fof_allocate();
    tree_build( pt );
    for ( i=0; i<npart; i++ ) {
        Head[i] = Tail[i] = i;
        Len[i] = 1;
        Next[i] = -1;
    }
    rhodm = header.Omega0 * 3 * SQR( header.HubbleParam ) / ( 8 * M_PI * All.G );
    if ( header.mass[pt] != 0 )
        mass = header.mass[pt];
    else{
        for ( i=offset, masstot=0; i<offset+npart; i++ )
            masstot += P[i].Mass;
        mass = masstot / npart;
    }
    LinkL = All.LinkLength * pow( mass / rhodm, 1.0/3 );
    sprintf( LogBuf, "critical density of dark matter: %g\n"
            "Comoving linking lenght %g", rhodm, LinkL );
    print_log( LogBuf );
    fof_find_groups();
    i = 0;
    i = Head[i];
    while( i>=0 ) {
        printf( "%li\n", i );
        i = Next[i];
    }
    tree_free();
    print_log( "fof ... done" );
    print_log( sep_str );
}
