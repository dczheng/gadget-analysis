#include "allvars.h"
#define SLICE_NUM 6


void generate_2D_img( char *fn_prefix, PLFLT **z, PLINT *nxy, PLFLT *xymm  ){
    char buf[50];
    PLFLT zmin, zmax;
    PLINT i,j;
    sprintf( buf, "%s.svg", fn_prefix );
    plsdev( "svg" );
    plsfnam( buf );
    sprintf( buf, "z=%.2f", redshift );
    plinit();
    plcol0( 2 );
    plenv( 0.1, 0.9, 0.1, 0.9, 1, -1 );
    pllab( "", "", buf );
    zmin = zmax = z[0][0];
    //zmin = zmax = 0;
    for ( i=0; i<nxy[0]; i++ )
        for ( j=0; j<nxy[1]; j++ ){
            if ( zmin==0  && z[i][j]!=0 ) zmin = zmax = z[i][j];
            zmin = ( z[i][j]<zmin ) ? z[i][j] : zmin;
            zmax = ( z[i][j]>zmax ) ? z[i][j] : zmax;
            /*
            zmin = ( z[i][j]<zmin && z[i][j] != 0 ) ? z[i][j] : zmin;
            zmax = ( z[i][j]>zmax ) ? z[i][j] : zmax;
            */
        }
    /*
    for ( i=0; i<nxy[0]; i++ )
        for ( j=0; j<nxy[1]; j++ ){
            if ( z[i][j] == 0 ) z[i][j] = zmin;
            z[i][j] = log10( z[i][j] );
        }
    zmin = log10( zmin );
    zmax = log10( zmax );
    */
    fprintf( stdout, "zmin=%e, zmax=%e\n", zmin, zmax );
    plimage( (PLFLT_MATRIX)z, nxy[0], nxy[1],
            0.0, 1.0, 0.0, 1.0,
            zmin, zmax,
            0.0, 1.0, 0.0, 1.0 );
    plend();
}

int compare_for_plot_baryon_density( const void *a, const void *b ) {
    float *aa, *bb;
    aa = ( float* )a;
    bb = ( float* )b;
    if ( aa[2] > bb[2] ) return 1;
    return -1;
}

void plot_scalar( int pt, enum iofields blk ){
    float *data, dz, r, *p;
    char fn_buf[50], buf[20];
    long i, N, j, z1, z2, k, test_num, ii, jj;
    PLINT nxy[2];
    PLFLT **rho, xymm[4];
    switch ( blk ) {
        case IO_U:
            if ( pt != 0 ) {
                get_dataset_name( blk, buf );
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            p = Particle[pt].u;
            break;
        case IO_RHO:
            if ( pt != 0 ) {
                get_dataset_name( blk, buf );
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            p = Particle[pt].rho;
            break;
        case IO_MASS:
            p = Particle[pt].m;
            break;
        case IO_POT:
            p = Particle[pt].m;
            break;
        case IO_ELEC:
            if ( pt != 0 ) {
                get_dataset_name( blk, buf );
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            p = Particle[pt].elec;
            break;
        case IO_VEL:
            p = ( float* ) malloc( sizeof(float) * Particle[pt].num );
            for ( i=0; i<Particle[pt].num; i++ ) {
                p[i] = sqrt( pow( Particle[pt].vel[i*3+0], 2 ) +
                             pow( Particle[pt].vel[i*3+1], 2 ) +
                             pow( Particle[pt].vel[i*3+2], 2 ) );
            }
            break;
        case IO_ACCEL:
            p = ( float* ) malloc( sizeof(float) * Particle[pt].num );
            for ( i=0; i<Particle[pt].num; i++ ) {
                p[i] = sqrt( pow( Particle[pt].accel[i*3+0], 2 ) +
                             pow( Particle[pt].accel[i*3+1], 2 ) +
                             pow( Particle[pt].accel[i*3+2], 2 ) );
            }
            break;
        case IO_MAG:
            if ( pt != 0 ) {
                get_dataset_name( blk, buf );
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            p = ( float* ) malloc( sizeof(float) * Particle[pt].num );
            for ( i=0; i<Particle[pt].num; i++ ) {
                p[i] = sqrt( pow( Particle[pt].mag[i*3+0], 2 ) +
                             pow( Particle[pt].mag[i*3+1], 2 ) +
                             pow( Particle[pt].mag[i*3+2], 2 ) );
            }
            break;
    }
    nxy[0] = pic_xsize;
    nxy[1] = pic_ysize;
    plAlloc2dGrid( &rho, nxy[0], nxy[1] );
    N = Particle[0].num;
    fputs( sep_str, stdout );
    get_dataset_name( blk, buf );
    fprintf( stdout, "plot partilcle %i field: \"%s\" ... \n", pt, buf );
    data = ( float* ) malloc( sizeof( float ) * N * 4 );
    for ( i=0; i<N; i++ ){
        data[ i*4+0 ] = Particle[pt].pos[ i*3+0 ];
        data[ i*4+1 ] = Particle[pt].pos[ i*3+1 ];
        data[ i*4+2 ] = Particle[pt].pos[ i*3+2 ];
        data[ i*4+3 ] = p[ i ];
    }
    test_num = 10;
    fputs( sep_str, stdout );
    for ( i=0; i<test_num; i++ ) {
        fprintf( stdout, "( %15.5f %15.5f %15.5f ):   %e\n",
                data[ i*4+0 ], data[ i*4+1 ], data[ i*4+2 ],
                data[ i*4+3 ] );
    }
    qsort( (void*)data, N, sizeof( float )*4,
                &compare_for_plot_baryon_density );
    fputs( "Sorted by z: \n", stdout );
    for ( i=0; i<test_num; i++ ) {
        fprintf( stdout, "( %15.5f %15.5f %15.5f ):   %e\n",
                data[ i*4+0 ], data[ i*4+1 ], data[ i*4+2 ],
                data[ i*4+3 ] );
    }
    fputs( sep_str, stdout );
    dz = header.BoxSize / slice_num;
    fprintf( stdout, "slice info: dz=%f\n", dz );
    for ( i=0; i<slice_num; i++ ) {
        for ( j=0; j<slice_index_num; j++ ) {
            if ( i == slice_index[j] ) {
                z1 = 0;
                while ( data[ z1*4+2 ] < i * dz ) z1++;
                z2 = z1;
                while ( data[ z2*4+2 ] < ( i+1 )*dz && z2<Particle[pt].num ) {
                    //fprintf( stdout, "%i: %f\n", z2, data[ z2*4+2 ] );
                    z2++;
                }
                z2--;
                fprintf( stdout, "slice_info: z1=%li, z2=%li ( %.4lf ~ %.4lf )\n",
                        z1, z2, i*dz, (i+1)*dz );
                for ( k=z1; k<z1+test_num; k++ ) {
                    fprintf( stdout, "( %15.5f %15.5f %15.5f ):   %e\n",
                            data[ k*4+0 ], data[ k*4+1 ], data[ k*4+2 ],
                            data[ k*4+3 ] );
                }
                for ( k=z2; k>z2-test_num; k-- ) {
                    fprintf( stdout, "( %15.5f %15.5f %15.5f ):   %e\n",
                            data[ k*4+0 ], data[ k*4+1 ], data[ k*4+2 ],
                            data[ k*4+3 ] );
                }
                for ( ii=0; ii< nxy[0]; ii++ )
                    for ( jj=0; jj<nxy[1]; jj++ )
                        rho[ii][jj] = 0;
                for ( k=z2; k>z1; k-- ) {
                    ii = ( PLINT )(data[ k*4+0 ] / ( header.BoxSize / nxy[0] ));
                    jj = ( PLINT )(data[ k*4+1 ] / ( header.BoxSize / nxy[1] ));
                    r = ( dz - (data[ k*4+2 ] - i*dz) ) / dz;
                    //r = ( float )( z2-k ) / ( z2-z1-1 ) ;
                    //fprintf( stdout, "%f\n", r );
                    rho[ii][jj] += ( PLFLT ) ( data[ k*4+3 ] * pow( r,1 ) );
                    //r=1;
                    //rho[ii][jj] += ( PLFLT ) ( data[ k*4+3 ] * pow( r,1 ) );
                }
                sprintf( fn_buf, "%s_%i_%i", Out_Picture_Prefix, ( int )( i*dz ), ( int )((i+1)*dz) );
                xymm[0] = i * dz;
                xymm[1] = (i+1) * dz;
                xymm[2] = i * dz;
                xymm[4] = (i+1) * dz;
                generate_2D_img( fn_buf, rho, nxy, xymm );
            }
        }
    }
    plFree2dGrid( rho, nxy[0], nxy[1] );
    fputs( sep_str, stdout );
    free( data );
    switch ( blk ) {
        case IO_MAG:
            if ( pt==0 ) free( p );
            break;
        case IO_VEL:
        case IO_ACCEL:
            free( p );
            break;
    }
}

int Compare_For_Plot_2D_Point( const void *a, const void *b ){
    return ( ( (float*)b )[2] < ( ( float*)a )[2] ) ? 1 : 0;
}

void Plot_2D_Point( int pt ) {
    PLFLT *x, *y, dz, z, nz, xmin, xmax, ymin, ymax;
    PLINT i, j, k, num, array_len, index, flag, s;
    PLINT slice_id[SLICE_NUM] = {
        10,
        20,
        30,
        40,
        50
    };
    char buf[50];
    float *pos;
    fputs( sep_str, stdout );
    fputs( "Plot 2D Point ...\n", stdout );
    pos = ( float* ) malloc( sizeof( float ) * Particle[pt].num * 3 );
    memcpy( pos, Particle[pt].pos, sizeof( float ) * Particle[pt].num *3 );
    /*
    for ( i=0; i<100; i++ ) {
        fprintf( stdout, "%f %f %f\n", pos[i*3+0], pos[i*3+1], pos[i*3+2] );
    }
    */
    fprintf( stdout, "\n" );
    qsort( (void*)pos, Particle[pt].num, sizeof( float )*3,
                &Compare_For_Plot_2D_Point );
    /*
    for ( i=0; i<100; i++ ) {
        fprintf( stdout, "%f %f %f\n", pos[i*3+0], pos[i*3+1], pos[i*3+2] );
    }
    */
    nz = 1000;
    dz = header.BoxSize / nz;
    z = dz;
    i = j = 0;
    index = 0;
    while( j<Particle[pt].num ) {
        while( pos[j*3+2]<z  &&
                j<Particle[pt].num ) j++;
        index ++;
        //fprintf( stdout, "i=%i j=%i num=%i\n", i, j, num  );
        flag = 0;
        for ( s=0; s<SLICE_NUM; s++ ) {
            if ( index == slice_id[s] ){
                flag = 1;
                break;
            }
        }
        if ( flag ) {
            num = j-i;
            x = ( PLFLT* ) malloc( sizeof( PLFLT ) * num );
            y = ( PLFLT* ) malloc( sizeof( PLFLT ) * num );
            for ( k=i; k<j; k++ ) {
                x[k-i] = pos[k*3+0];
                y[k-i] = pos[k*3+1];
            }
            plsdev( "pngcairo" );
            sprintf( buf, "%s_%i_%i.png", Out_Picture_Prefix,
                        (int)(z-dz), (int)z);
            plsfnam( buf );
            plinit();
            xmin = 0;
            xmax = header.BoxSize;
            ymin = 0;
            ymax = header.BoxSize;
            plenv( xmin, xmax, ymin, ymax, 0,0 );
            plpoin( num, x, y, 1 );
            plend();
            free( x );
            free( y );
        }
        i = j;
        z += dz;
    }
    free( pos );
    fputs( sep_str, stdout );
}
