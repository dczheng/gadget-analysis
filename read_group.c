#include "allvars.h"

void read_group() {
    DIR *dir;
    FILE *fd;
    char group_file[ FILENAME_MAX ];
    struct dirent *ptr;
    struct group_struct *group_local;
    int Ngroups, Ntask, Nids, malloc_flag;
    int group_offset, *len, i, j, *offset, *lentype;
    long TotNids;
    float *mass, *cm, *vel, *veldisp, *tensor, *rmax, *vmax, *pos;
    float *angmom;
    fputs( sep_str, stdout );
    fputs( "read group ...\n", stdout );
    if ( !( dir=opendir( All.GroupDir ) ) ){
        fprintf( stderr, "Failed to open GroupDir %s\n", All.GroupDir );
        endrun( 4 );
    }
    malloc_flag = 1;
    group_offset = 0;
    while (  ptr=readdir( dir ) ) {
        if ( ( strcmp( ptr->d_name, "." ) == 0 ) ||
             ( strcmp( ptr->d_name, ".." ) == 0 ) )
            continue;
        sprintf( group_file, "%s%s", All.GroupDir, ptr->d_name );
#ifdef GROUP_ZDEBUG
        fprintf( stdout, "%s\n", ptr->d_name );
#endif
        if ( !( fd = fopen( group_file, "r" ) ) ) {
            fprintf( stderr, "Failed to open group file :%s\n", group_file );
            endrun( 5 );
        }
        fread( &Ngroups, sizeof( int ), 1, fd );
        fread( &TotNgroups, sizeof( int ), 1, fd );
        fread( &Nids, sizeof( int ), 1, fd );
        fread( &TotNids, sizeof( long ), 1, fd );
        fread( &Ntask, sizeof( int ), 1, fd );
        if ( malloc_flag ) {
            group = ( struct group_struct * ) malloc(
                    sizeof( struct group_struct ) *
                    TotNgroups );
            malloc_flag = 0;
        }
        if ( Ngroups ){
#ifdef GROUP_ZDEBUG
            fprintf( stdout, "Ngroups = %i, TotNgroups = %i,"
                    " Nids = %i, TotNids = %lli, Ntask = %i\n",
                    Ngroups, TotNgroups, Nids, TotNids, Ntask );
#endif

            len = ( int * ) malloc( sizeof( int ) * Ngroups );
            if ( ! fread( len, sizeof( int ), Ngroups, fd ) ) {
                fprintf( stderr, "Failed to read Len array!\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ) {
#ifdef GROUP_ZDEBUG
                fprintf( stdout, "len[%i] = %i\n", i, len[i] );
#endif
                group[i+group_offset].Len = len[i];
            }
            free( len );

            offset = ( int * ) malloc( sizeof( int ) * Ngroups );
            if ( !fread( offset, sizeof( int ), Ngroups, fd ) ){
                fprintf( stderr, "Failed to read Offset array!\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ) {
#ifdef GROUP_ZDEBUG
                fprintf( stdout, "offset[%i] = %i\n", i, offset[i] );
#endif
                group[i+group_offset].Offset = offset[i];
            }
            free( offset );

            mass = ( float* ) malloc( sizeof( float ) * Ngroups );
            if ( !fread( mass, sizeof( float ), Ngroups, fd ) ) {
                fprintf( stderr, "Failed to read Mass array!\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ) {
#ifdef GROUP_ZDEBUG
                fprintf( stdout, "mass[%i] = %f\n", i, mass[i] );
#endif
                group[i+group_offset].Mass = mass[i];
            }
            free( mass );

            cm = ( float* ) malloc( sizeof( float ) * Ngroups * 3 );
            if ( !fread( cm, sizeof( float )*3, Ngroups, fd ) ){
                fprintf( stderr, "Failed to read CM array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ){
                for ( j=0; j<3; j++ ) {
#ifdef GROUP_ZDEBUG
                    fprintf( stdout, "cm[%i][%i]=%-15.7f ", i, j, cm[ i*3+j ] );
#endif
                    group[ i+group_offset ].CM[j] = cm[ i*3+j ];
                }
#ifdef GROUP_ZDEBUG
                fprintf( stdout, "\n" );
#endif
            }
            free( cm );

            vel = ( float* ) malloc( sizeof( float ) * Ngroups * 3 );
            if ( !fread( vel, sizeof( float )*3, Ngroups, fd ) ) {
                fprintf( stderr, "Failed to read Vel array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ){
                for ( j=0; j<3; j++ ) {
#ifdef GROUP_ZDEBUG
                    fprintf( stdout, "vel[%i][%i]=%-15.7f ", i, j, vel[ i*3+j ] );
#endif
                    group[ i+group_offset ].Vel[j] = vel[ i*3+j ];
                }
#ifdef GROUP_ZDEBUG
                fprintf( stdout, "\n" );
#endif
            }
            free( vel );

#ifndef FOF_EXTENDED_PROPERTIES
            lentype = ( int * ) malloc( sizeof( int ) * Ngroups * 6 );
            if ( !fread( lentype, sizeof( int )*6, Ngroups, fd ) ) {
                fprintf( stderr, "Failed to read LenType array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ) {
                for ( j=0; j<6; j++ ) {
#ifdef GROUP_ZDEBUG
                    fprintf( stdout, "lentype[%i][%i]=%i ", i, j, lentype[ i*6+j ] );
#endif
                    group[ i+group_offset ].LenType[j] = len[ i*6+j ];
                }
#ifdef GROUP_ZDEBUG
                fprintf( stdout, "\n" );
#endif
            }
            free( lentype );
#endif
#ifdef FOF_EXTENDED_PROPERTIES
            veldisp = ( float* ) malloc( sizeof( float ) * Ngroups );
            if ( !fread( veldisp, sizeof( float ), Ngroups ) ) {
                fprintf( stderr, "Failed to read VelDisp array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ) {
#ifdef GROUP_ZDEBUG
                fprintf( stdout, "veldisp[i]=%f\n", i, veldisp[i] );
#endif
                group[ i+group_offset ].VelDisp = veldisp[i];
            }
            free( veldisp );

            tensor = ( float* ) malloc( sizeof( float ) * Ngroups * 9 );
            if( !fread( tensor, sizeof( float ) * 9, Ngroups, fd ) ) {
                fprintf( stderr, "Failed to read ToI array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ) {
                for( j=0; j<9; j++ ) {
                    group[ i+group_offset ].ToI[j] = tensor[ i*9+j ];
                }
            }
            free( tensor );

            rmax = ( float* ) malloc( sizeof( float ) * Ngroups );
            if( !fread( rmax, sizeof( float ), Ngroups, fd) ) {
                fprintf( stderr, "Failed to read Rmax array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ){
                group[ i+group_offset ].Rmax = rmax[i];
            }
            free( rmax );

            vmax = ( float* ) malloc( sizeof( float ) * Ngroups );
            if( !fread( vmax, sizeof( float ), Ngroups, fd) ) {
                fprintf( stderr, "Failed to read Vmax array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ){
                group[ i+group_offset ].Vmax = vmax[i];
            }
            free( rmax );

            pos = ( float* ) malloc( sizeof( float ) * Ngroups * 3 );
            if( !fread( pos, sizeof( float )*3, Ngroups, fd) ) {
                fprintf( stderr, "Failed to read Pos array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ){
                for ( j=0; j<3; j++ ) {
                    group[ i+group_offset ].Pos[j] = pos[ i*3+j ];
                }
            }
            free( pos );

            angmom = ( float* ) malloc( sizeof( float ) * Ngroups * 3 );
            if( !fread( angmom, sizeof( float )*3, Ngroups, fd) ) {
                fprintf( stderr, "Failed to read AngMom array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ){
                for ( j=0; j<3; j++ ) {
                    group[ i+group_offset ].AngMom[j] = angmom[ i*3+j ];
                }
            }
            free( angmom );
#endif
#ifdef SFR
            mass = ( float* ) malloc( sizeof( float ) * Ngroups );
            if( !fread( mass, sizeof( float ), Ngroups, fd ) ) {
                fprintf( stderr, "Failed to read Sfr array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ) {
#ifdef GROUP_ZDEBUG
                fprintf( stdout, "sfr_mass[%i]=%f\n", i, mass[i] );
#endif
                group[ i+group_offset ].Sfr = mass[i];
            }
            free( mass );
#endif
#ifdef BLACK_HOLES
            mass = ( float* ) malloc( sizeof( float ) * Ngroups );
            if( !fread( mass, sizeof( float ), Ngroups, fd ) ) {
                fprintf( stderr, "Failed to read BH_Mass array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ) {
#ifdef GROUP_ZDEBUG
                fprintf( stdout, "bh_mass[%i]=%f\n", i, mass[i] );
#endif
                group[ i+group_offset ].BH_Mass = mass[i];
            }
            free( mass );

            mass = ( float* ) malloc( sizeof( float ) * Ngroups );
            if ( !fread( mass, sizeof( float ), Ngroups, fd ) ) {
                fprintf( stderr, "Failed to read BH_Mdot array.\n" );
                endrun( 6 );
            }
            for ( i=0; i<Ngroups; i++ ) {
#ifdef GROUP_ZDEBUG
                fprintf( stdout, "bh_mdot_mass[%i]\n", i, mass[i] );
#endif
                group[ i+group_offset ].BH_Mdot = mass[i];
            }
            free( mass );
#endif
            group_offset += Ngroups;
        }
        fclose( fd );
    }
    closedir( dir );
    fputs( sep_str, stdout );
}

void free_group() {
    free( group );
}
