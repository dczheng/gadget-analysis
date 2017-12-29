#include "allvars.h"

long npart, offset;
double LinkL, rhodm;

void fof_allocate() {
    size_t bytes;

    if ( !(fof_info = malloc( bytes = npart * sizeof( struct fof_info_struct ) ) ) )
        printf( "failed to allocate memory for `fof_info` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `fof_info` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(fof_Next = malloc( bytes = npart * sizeof( long ) ) ) )
        printf( "failed to allocate memory for `fof_Next` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `fof_Next` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    if ( !(Ngblist = malloc( bytes = npart * sizeof( long ) ) ) )
        printf( "failed to allocate memory for `Ngblist` ( %g MB )\n", bytes / 1024.0 / 1024.0 );
    sprintf( LogBuf, "allocate memory for `Ngblist` ( %g MB )", bytes / 1024.0 / 1024.0 );
    print_log( LogBuf );

    print_log( sep_str );
}

void fof_free() {
    free( fof_info );
    free( fof_Next );
    free( Ngblist );
}

int fof_compare_len( const void *a, const void *b ) {
    return ( ( ( ( struct fof_info_struct *)a )->Len < ( ( struct fof_info_struct * )b )->Len ) ? 1 : -1 );
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
            if ( fof_info[ii].Head != fof_info[j].Head ) {
                if ( fof_info[fof_info[ii].Head].Len > fof_info[fof_info[j].Head].Len ) {
                    p = ii;
                    s = j;
                }
                else{
                    p = j;
                    s = ii;
                }
                fof_Next[ fof_info[fof_info[p].Head].Tail ] = fof_info[s].Head;
                fof_info[ fof_info[p].Head ].Tail = fof_info[ fof_info[s].Head ].Tail;
                fof_info[ fof_info[p].Head ].Len += fof_info[ fof_info[s].Head ].Len;
                ss = fof_info[s].Head;
                do
                    fof_info[ss].Head = fof_info[p].Head;
                while( (ss=fof_Next[ss]) >=0 );
            }
        }
    }
    qsort( fof_info, npart, sizeof( struct fof_info_struct ), fof_compare_len );
    print_log( "fof find groups ... done" );
}

void fof_compute_group_properties() {
    int i;
    print_log( "fof compute groups properties ... " );
    for ( i=0; i<npart; i++ ) {
        if ( fof_info[i].Len > All.FofMinLen )
            printf( "%li\n", fof_info[i].Len );
    }
    print_log( "fof compute groups properties ... done " );
}

void fof_test() {
    int i, index;
    long p, j;
    FILE *fd;
    index = 0;
    fd = fopen( "fof.txt", "w" );
    p = fof_info[index].Head;
    for ( i=0; i<fof_info[index].Len; i++ ){
        j = p + offset;
        fprintf( fd, "%g %g %g\n",
                P[j].Pos[0],
                P[j].Pos[1],
                P[j].Pos[2] );
        p = fof_Next[p];
    }
    fclose( fd );
}

void fof_save() {
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
        fof_info[i].Head = fof_info[i].Tail = i;
        fof_info[i].Len = 1;
        fof_Next[i] = -1;
    }
    rhodm = (All.Omega0-All.OmegaBaryon) * 3 * SQR( All.Hubble ) / ( 8 * M_PI * All.G );
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
    fof_compute_group_properties();
    fof_test();
    tree_free();
    print_log( "fof ... done" );
    print_log( sep_str );
}
