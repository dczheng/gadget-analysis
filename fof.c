#include "allvars.h"

long npart, offset, num;
MyIDType *Head, *Len, *Next, *Tail;
int Ngroups;
void fof_allocate() {
    size_t bytes;
    if ( !(Head = malloc( bytes = num * sizeof( MyIDType ) ) ) )
        printf( "failed to allocate memory for `Head` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `Head` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(Tail = malloc( bytes = num * sizeof( MyIDType ) ) ) )
        printf( "failed to allocate memory for `Tail` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `Tail` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(Len = malloc( bytes = num * sizeof( MyIDType ) ) ) )
        printf( "failed to allocate memory for `Len` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `Len` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(Next = malloc( bytes = num * sizeof( MyIDType ) ) ) )
        printf( "failed to allocate memory for `Next` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `Next` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(Ngblist = malloc( bytes = num * sizeof( int ) ) ) )
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
    print_log( "start fof find groups ..." );
    print_log( "fof find groups ... done" );
}

void fof( int pt ) {
    long i;
    print_log( "start fof ..." );
    num = get_particle_num( pt );
    offset = get_particle_offset( pt );
    for ( i=offset, num=0; i<offset+npart; i++ )
        if ( P[i].Flag ) num++;
    fof_allocate();
    tree_build( pt );
    fof_find_groups();
    tree_free();
    print_log( "fof ... done" );
    print_log( sep_str );
}
