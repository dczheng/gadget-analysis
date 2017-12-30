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

void fof_find_groups( int pt ) {
    int ngbnum, k, ngbmax;
    double pos[3];
    long i, ii, p, s, j, ss;
    //int flag=0;
    ngbmax = 0;
    print_log( "start fof find groups ..." );
    for ( i=offset; i<offset+npart; i++ ) {
        ii = i-offset;
        if ( !P[i].Flag ) continue;
        for( j=0; j<3; j++ )
            pos[j] = P[i].Pos[j];
        ngbnum = ngb_fof( pos, LinkL, pt );
        //printf( "%i\n", ngbnum );
        //if ( flag++ > 30 ) endrun( 20171230 );
        if ( ngbnum > ngbmax ) ngbmax = ngbnum;
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
                fof_info[ fof_info[s].Head ].Len = 0;
                ss = fof_info[s].Head;
                do
                    fof_info[ss].Head = fof_info[p].Head;
                while( (ss=fof_Next[ss]) >= 0 );
            }
        }
    }
    qsort( fof_info, npart, sizeof( struct fof_info_struct ), fof_compare_len );
    sprintf( LogBuf, "the maximum number of ngb: %i", ngbmax );
    print_log( LogBuf );
    print_log( "fof find groups ... done" );
}

void fof_compute_group_properties( int pt ) {
    int i, num, j, k;
    long ii, p;
    double mass;
    print_log( "fof compute groups properties ... " );
    Ngroups = 0;
    for ( i=0; i<npart; i++ ) {
        if ( fof_info[i].Len < All.FofMinLen ) break;
        p = fof_info[i].Head;
        fof_info[i].mass = 0;
        for ( k=0; k<3; k++ ){
            fof_info[i].cm[k] = 0;
            fof_info[i].vel[k] = 0;
        }
        fof_info[i].vr200 = 0;
        for ( j=0; j<fof_info[i].Len; j++ ) {
            ii = p + offset;
            p = fof_Next[p];
            mass = ( header.mass[pt] != 0 ) ? header.mass[pt] : P[ii].Mass;
            fof_info[i].mass += mass;
            for ( k=0; k<3; k++ ){
                fof_info[i].cm[k] += mass * P[ii].Pos[k];
                fof_info[i].vel[k] += mass * P[ii].Vel[k];
            }
        }
        for ( k=0; k<3; k++ ){
            fof_info[i].cm[k] /= fof_info[i].mass;
            fof_info[i].vel[k] /= fof_info[i].mass;
        }
        fof_info[i].vr200 = pow( fof_info[i].mass /
            ( All.CriticalDensity * 200 * 4.0 / 3.0 * PI ), 1.0/3.0 );
    }
    Ngroups = i;
    sprintf( LogBuf, "The number of groups with at least %i particles: %li", All.FofMinLen, Ngroups );
    print_log( LogBuf );
    sprintf( LogBuf, "Largest group has %li particles", fof_info[0].Len );
    print_log( LogBuf );
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
    hid_t hdf5_file, hdf5_dataset, hdf5_dataspace;
    herr_t herr;
    sprintf( All.FofFileName, "%s.hdf5", All.FofFileName );
    hdf5_file = H5Fcreate( All.FofFileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    H5Fclose( hdf5_file );
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
    rhodm = (All.Omega0-All.OmegaBaryon) * 3 * SQR( All.Hubble ) / ( 8 * PI * All.G );
    if ( header.mass[pt] != 0 )
        mass = header.mass[pt];
    else{
        for ( i=offset, masstot=0; i<offset+npart; i++ )
            masstot += P[i].Mass;
        mass = masstot / npart;
    }
    LinkL = All.LinkLength * pow( mass / rhodm, 1.0/3 );
    sprintf( LogBuf, "critical density of dark matter: %g\n"
            "comoving linking lenght %g", rhodm, LinkL );
    print_log( LogBuf );
    fof_find_groups( pt );
    fof_compute_group_properties( pt );
    fof_test();
    tree_free();
    print_log( "fof ... done" );
    print_log( sep_str );
}
