#include "allvars.h"

void generate_2D_img( char *fn_prefix, char *x_buf, char *y_buf, PLFLT **z, PLINT *nxy,
        int log_flag, PLFLT *xymm, PLFLT **u, PLFLT **v, PLcGrid2 *cgrid2, PLINT *anxy, int flag_a ){
    char buf[50];
    PLFLT zmin, zmax, colorbar_width, colorbar_heigh;
    PLINT i,j;
    PLFLT arr_x[6] = { -0.5, 0.5, 0.3, 0.5, 0.3, 0.5 };
    PLFLT arr_y[6] = { 0.0, 0.0, 0.2, 0.0, -0.2, 0.0 };
    PLINT arr_n = 6;
    sprintf( buf, "%s.png", fn_prefix );
    plsdev( "png" );
    plsfnam( buf );
    sprintf( buf, "z=%.2f", redshift );
    plspage( 100, 100, 1024, 1024, 0, 0 );
    plinit();
    plcol0( 15 );
    /*
    plvpor( 0.0, 1.0, 0.0, 1.0 );
    plvpor( 0,0, (PLFLT)header.BoxSize/1000,
            0.0, (PLFLT)header.Boxsize/1000 );
    plbox( "n", 0, 0, "n", 0, 0 );
           */
    plenv( xymm[0]/1000, xymm[1]/1000 ,
           xymm[2]/1000, xymm[3]/1000 ,
           1, 0 );
    pllab( x_buf, y_buf, buf );
    //zmin = zmax = z[0][0];
    zmin = zmax = 0;
    for ( i=0; i<nxy[0]; i++ )
        for ( j=0; j<nxy[1]; j++ ){
            //zmin = ( z[i][j]<zmin ) ? z[i][j] : zmin;
            //zmax = ( z[i][j]>zmax ) ? z[i][j] : zmax;
            //fprintf( stdout, "%e\n", z[i][j] );
            if ( zmin==0  && z[i][j]!=0 ) zmin = zmax = z[i][j];
            zmin = ( z[i][j]<zmin && z[i][j] != 0 ) ? z[i][j] : zmin;
            zmax = ( z[i][j]>zmax ) ? z[i][j] : zmax;
        }
    zmin = zmin / 10.0;
    if ( log_flag ) {
        for ( i=0; i<nxy[0]; i++ )
            for ( j=0; j<nxy[1]; j++ ){
                if ( z[i][j] == 0 ) z[i][j] = zmin;
                z[i][j] = log10( z[i][j] );
            }
        zmin = log10( zmin );
        zmax = log10( zmax );
    }
    fprintf( stdout, "zmin=%e, zmax=%e\n", zmin, zmax );
    plspal1( "./plot.pal", 1 );
    plimage( (PLFLT_MATRIX)z, nxy[0], nxy[1],
            xymm[0]/1000, xymm[1]/1000,
            xymm[2]/1000, xymm[3]/1000,
            zmin, zmax,
            xymm[0]/1000, xymm[1]/1000,
            xymm[2]/1000, xymm[3]/1000 );
    //set colorbar
    PLINT n_axis = 1;
    PLCHAR_VECTOR axis_optsl[] = { "vstml", };
    PLCHAR_VECTOR axis_opts[] = { "vstm", };
    PLFLT axis_ticks[1] = { 0.0, };
    PLINT axis_subticks[1] = { 0, };
    PLINT num_values[1] = { 15 };
    PLFLT *values[1];
    values[0] = ( PLFLT* ) malloc( num_values[0] * sizeof( PLFLT ) );
    for ( i=0; i<num_values[0]; i++ ) {
        values[0][i] = i * ( zmax-zmin ) / num_values[0] + zmin;
    }
    plcol0( 2 );
    plschr( 0, 0.5 );
    if ( log_flag )
        plcolorbar( &colorbar_width, &colorbar_heigh,
            PL_COLORBAR_IMAGE, 0,
            0.02, 0, 0.03, 0.9, 0, 1, 1, 0.0, 0.0, 0.0, 0.0,
            0, NULL, NULL,
            n_axis, axis_optsl,
            axis_ticks, axis_subticks,
            num_values, (PLFLT_MATRIX) values );
    else
        plcolorbar( &colorbar_width, &colorbar_heigh,
            PL_COLORBAR_IMAGE, 0,
            0.02, 0, 0.03, 0.9, 0, 1, 1, 0.0, 0.0, 0.0, 0.0,
            0, NULL, NULL,
            n_axis, axis_opts,
            axis_ticks, axis_subticks,
            num_values, (PLFLT_MATRIX) values );
    free( values[0] );
    //set colorbar
    if ( flag_a ) {
    plcol0( 15 );
if ( this_task == 0 )
    fputs( "plot arrow ...\n", stdout );
    plsvect( arr_x, arr_y, arr_n, 0 );
    plvect( (PLFLT_MATRIX)u, ( PLFLT_MATRIX )v, anxy[0], anxy[1], 0.0, pltr2, (void*)cgrid2 );
    }
    plend();
}

int compare_for_plot_slice( const void *a, const void *b ) {
    float *aa, *bb;
    aa = ( float* )a;
    bb = ( float* )b;
    if ( aa[2] > bb[2] ) return 1;
    return -1;
}

double kernel( double r, double h ) {
    double v;
    v = r/h;
    if ( ( 0<=v ) && ( v<=0.5 ) )
        return 8.0 / M_PI / pow( h, 3.0 ) * ( 1 - 6*pow( v, 2 ) + 6*pow( v, 3 ) ); //* r * r;
    if( ( 0.5<v ) && ( v<= 1) )
        return 8.0 / M_PI / pow( h, 3.0 ) * 2 * pow( 1-v, 3 ); // * r * r;
    /*
    if ( ( 0<=v ) && ( v<=1 ) )
        return 3.0 / 4.0 / M_PI / pow( h, 3 ) * (10.0/3.0 - 7 * pow( v, 2 ) + 4 * pow( v, 3 ));
    if( ( 1<v ) && ( v<= 2) )
        return 3.0 / 4.0 / M_PI / pow( h, 3 ) * ( pow( 2-v, 2 ) * ( 5-4*v ) / 3.0 );
    */
    return 0;
}

void plot_slice( int pt, enum iofields blk ){
    float *data, dz, r, *p, r1, r2;
    char fn_buf[50], buf[50], x_buf[100], y_buf[100], log_flag, mean_flag;
    long i, N, j, z1, z2, k, test_num, ii, jj, index, i1, i2, index1;
    long index_i, index_j, index_k, flag, index_i1, index_j1, index_k1, kN;
    double g, h, rr1, rr2, rr3, rr4, tmp, *kmatrix, g2, ag, max_mag;
    PLINT nxy[2], anxy[2];
    PLFLT **v3, xymm[4], **u, **v;
    PLcGrid2 cgrid2;
    double *v1, *v2;;
    xymm[0] = (PLFLT) slice_corner1[0];
    xymm[1] = (PLFLT) slice_corner2[0];
    xymm[2] = (PLFLT) slice_corner1[1];
    xymm[3] = (PLFLT) slice_corner2[1];
    nxy[0] = pic_xsize;
    nxy[1] = pic_ysize;
    anxy[0] = arrow_x;
    anxy[1] = arrow_y;
    cgrid2.nx = anxy[0];
    cgrid2.ny = anxy[1];
    switch ( blk ) {
        case IO_U:
            get_dataset_name( blk, buf );
            if ( pt != 0 ) {
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            sprintf( x_buf, "log(%s)", buf );
            log_flag = 1;
            mean_flag = 0;
            p = Particle[pt].u;
            break;
        case IO_RHO:
            get_dataset_name( blk, buf );
            if ( pt != 0 ) {
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            sprintf( x_buf, "log(%s)(10^10Msun/Kpc^2)", buf );
            log_flag = 1;
            mean_flag = 1;
            p = Particle[pt].rho;
            plAlloc2dGrid( &u, anxy[0], anxy[1] );
            plAlloc2dGrid( &v, anxy[0], anxy[1] );
            plAlloc2dGrid( &cgrid2.xg, anxy[0], anxy[1] );
            plAlloc2dGrid( &cgrid2.yg, anxy[0], anxy[1] );
            break;
        case IO_MASS:
            get_dataset_name( blk, buf );
            p = Particle[pt].m;
            sprintf( x_buf, "log(%s)(10^10Msun/Kpc^2)", buf );
            plAlloc2dGrid( &u, anxy[0], anxy[1] );
            plAlloc2dGrid( &v, anxy[0], anxy[1] );
            plAlloc2dGrid( &cgrid2.xg, anxy[0], anxy[1] );
            plAlloc2dGrid( &cgrid2.yg, anxy[0], anxy[1] );
            log_flag = 1;
            mean_flag = 1;
            break;
        case IO_POT:
            get_dataset_name( blk, buf );
            if ( pt != 0 ) {
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            sprintf( x_buf, "log(%s)", buf );
            log_flag = 1;
            mean_flag = 0;
            p = Particle[pt].pot;
            break;
        case IO_MN:
            get_dataset_name( blk, buf );
            if ( pt != 0 ) {
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            sprintf( x_buf, "%s", buf );
            log_flag = 0;
            mean_flag = 0;
            p = Particle[pt].mn;
            break;
        case IO_J:
            get_dataset_name( blk, buf );
            if ( pt != 0 ) {
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            sprintf( x_buf, "log(%s)", buf );
            log_flag = 1;
            mean_flag = 0;
            p = Particle[pt].j;
            break;
        case IO_ELEC:
            get_dataset_name( blk, buf );
            if ( pt != 0 ) {
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            sprintf( x_buf, "%s", buf );
            log_flag = 0;
            mean_flag = 0;
            p = Particle[pt].elec;
            break;
        case IO_VEL:
            get_dataset_name( blk, buf );
            p = ( float* ) malloc( sizeof(float) * Particle[pt].num );
            for ( i=0; i<Particle[pt].num; i++ ) {
                p[i] = sqrt( pow( Particle[pt].vel[i*3+0], 2 ) +
                             pow( Particle[pt].vel[i*3+1], 2 ) +
                             pow( Particle[pt].vel[i*3+2], 2 ) );
            }
            sprintf( x_buf, "%s(km/s^2)", buf );
            plAlloc2dGrid( &u, anxy[0], anxy[1] );
            plAlloc2dGrid( &v, anxy[0], anxy[1] );
            plAlloc2dGrid( &cgrid2.xg, anxy[0], anxy[1] );
            plAlloc2dGrid( &cgrid2.yg, anxy[0], anxy[1] );
            log_flag = 0;
            mean_flag = 0;
            break;
        case IO_ACCEL:
            get_dataset_name( blk, buf );
            p = ( float* ) malloc( sizeof(float) * Particle[pt].num );
            for ( i=0; i<Particle[pt].num; i++ ) {
                p[i] = sqrt( pow( Particle[pt].accel[i*3+0], 2 ) +
                             pow( Particle[pt].accel[i*3+1], 2 ) +
                             pow( Particle[pt].accel[i*3+2], 2 ) );
            }
            log_flag = 0;
            mean_flag = 0;
            sprintf( x_buf, "%s", buf );
            break;
        case IO_MAG:
            get_dataset_name( blk, buf );
            if ( pt != 0 ) {
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            p = ( float* ) malloc( sizeof(float) * Particle[pt].num );
            for ( i=0; i<Particle[pt].num; i++ ) {
                p[i] = sqrt( pow( Particle[pt].mag[i*3+0], 2 ) +
                             pow( Particle[pt].mag[i*3+1], 2 ) +
                             pow( Particle[pt].mag[i*3+2], 2 ) ) / scalar_unit;
            }
            plAlloc2dGrid( &u, anxy[0], anxy[1] );
            plAlloc2dGrid( &v, anxy[0], anxy[1] );
            plAlloc2dGrid( &cgrid2.xg, anxy[0], anxy[1] );
            plAlloc2dGrid( &cgrid2.yg, anxy[0], anxy[1] );
            sprintf( x_buf, "log(%s)(uG/Kpc^2)", buf );
            log_flag = 1;
            mean_flag = 1;
            break;
    }
    plAlloc2dGrid( &v3, nxy[0], nxy[1] );
    v1 = ( double * ) malloc( sizeof( double ) * nxy[0] * nxy[1] );
if ( this_task == 0 )
        v2 = ( double * ) malloc( sizeof( double ) * nxy[0] * nxy[1] );
    //N = Particle[pt].num;
if ( this_task == 0 ){
    fputs( sep_str, stdout );
    get_dataset_name( blk, buf );
    fprintf( stdout, "plot partilcle %i field: \"%s\" ... \n", pt, buf );
}
    N = 0;
    max_mag = -1;
    for ( i=0; i<Particle[pt].num; i++ ) {
        if ( Particle[pt].pos[i*3+0] >= slice_corner1[0] &&
             Particle[pt].pos[i*3+0] <= slice_corner2[0] &&
             Particle[pt].pos[i*3+1] >= slice_corner1[1] &&
             Particle[pt].pos[i*3+1] <= slice_corner2[1] ){
            N++;
            switch ( blk ) {
                case IO_MAG:
                    if ( p[i] > max_mag ) max_mag = p[i];
                    break;
            }
        }
    }
if ( this_task == 0 )
    fprintf( stdout, "total point number: %li\n", N );
    data = ( float* ) malloc( sizeof( float ) * N * 6 );
    index = 0;
    for ( i=0; i<Particle[pt].num; i++ ){
        if ( Particle[pt].pos[i*3+0] >= slice_corner1[0] &&
             Particle[pt].pos[i*3+0] <= slice_corner2[0] &&
             Particle[pt].pos[i*3+1] >= slice_corner1[1] &&
             Particle[pt].pos[i*3+1] <= slice_corner2[1] ) {
                data[ index*6+0 ] = Particle[pt].pos[ i*3+0 ];
                data[ index*6+1 ] = Particle[pt].pos[ i*3+1 ];
                data[ index*6+2 ] = Particle[pt].pos[ i*3+2 ];
                data[ index*6+3 ] = p[ i ];
                switch ( blk ) {
                    case IO_MAG:
                        data[ index*6+4 ] = Particle[pt].mag[ i*3+0 ] / scalar_unit / max_mag;
                        /*
                            sqrt( pow( Particle[pt].mag[ i*3+0 ], 2 ) +
                                  pow( Particle[pt].mag[ i*3+1 ], 2 ) +
                                  pow( Particle[pt].mag[ i*3+2 ], 2 ) );
                                  */
                        data[ index*6+5 ] = Particle[pt].mag[ i*3+1 ] / scalar_unit / max_mag;
                        /*
                            sqrt( pow( Particle[pt].mag[ i*3+0 ], 2 ) +
                                  pow( Particle[pt].mag[ i*3+1 ], 2 ) +
                                  pow( Particle[pt].mag[ i*3+2 ], 2 ) );
                                  */
                        break;
                    case IO_MASS:
                    case IO_RHO:
                    case IO_VEL:
                        data[ index*6+4 ] = Particle[pt].vel[ i*3+0 ] / scalar_unit;
                        data[ index*6+4 ] = log10( data[ index*6+4 ] );
                        /*
                            sqrt( pow( Particle[pt].vel[ i*3+0 ], 2 ) +
                                  pow( Particle[pt].vel[ i*3+1 ], 2 ) +
                                  pow( Particle[pt].vel[ i*3+2 ], 2 ) );
                                  */
                        data[ index*6+5 ] = Particle[pt].vel[ i*3+1 ] / scalar_unit;
                        data[ index*6+5 ] = log10( data[ index*6+5 ] );
                        /*
                            sqrt( pow( Particle[pt].vel[ i*3+0 ], 2 ) +
                                  pow( Particle[pt].vel[ i*3+1 ], 2 ) +
                                  pow( Particle[pt].vel[ i*3+2 ], 2 ) );
                                  */
                        break;
                }
                index++;
        }
    }
    test_num = 10;
if ( this_task == 0 ){
        fputs( sep_str, stdout );
        for ( i=0; i<test_num; i++ ) {
            fprintf( stdout, "( %15.5f %15.5f %15.5f ):   %e\n",
                    data[ i*6+0 ], data[ i*6+1 ], data[ i*6+2 ],
                    data[ i*6+3 ] );
        }
}
    if ( slice_num > 1 ) {
        qsort( (void*)data, N, sizeof( float )*6,
                &compare_for_plot_slice );
        if ( this_task == 0 ){
            fputs( "Sorted by z: \n", stdout );
            for ( i=0; i<test_num; i++ ) {
                fprintf( stdout, "( %15.5f %15.5f %15.5f ):   %e\n",
                        data[ i*6+0 ], data[ i*6+1 ], data[ i*6+2 ],
                        data[ i*6+3 ] );
            }
        fputs( sep_str, stdout );
        }
    }
    dz = header.BoxSize / slice_num;
if ( this_task == 0 )
        fprintf( stdout, "slice info: dz=%f\n", dz );
    g = ( slice_corner2[0] - slice_corner1[0] ) / nxy[0];
    ag = ( slice_corner2[0] - slice_corner1[0] ) / anxy[0];
    h = 2.5;
    switch( blk ) {
        case IO_RHO:
            for ( i=0; i<N; i++ )
                data[ i*6+3 ] *= g * g * h;
            break;
    }
if ( this_task == 0 )
    fprintf( stdout, "g=%f, h=%f\n", g, h );
    if ( proj_mode == 1 ) {
if ( this_task == 0 )
    fputs( "initialize kernel matrix ...\n", stdout );
        kN = 200;
        kmatrix = ( double* ) malloc( kN * kN * sizeof( double ) );
        memset( kmatrix, 0,  kN * kN * sizeof(double) );
        for ( i=0; i<kN; i++ )
            for ( j=0; j<kN; j++ )
                for ( k=0; k<kN; k++ ) {
                    r = sqrt( pow( (i-kN/2)*h/kN, 2 ) +
                              pow( (j-kN/2)*h/kN, 2 ) +
                              pow( (k-kN/2)*h/kN, 2 ) );
                    kmatrix[ i*kN+j ] += kernel( r, h ) * pow( h/kN, 3 );
            }
    }
/*
if ( this_task == 0 ) {
    FILE *fd;
    fputs( "save kernel data to kernel.dat\n", stdout );
    fd = fopen( "kernel.dat", "w" );
    for ( i=0; i<kN; i++ ) {
        for ( j=0; j<kN; j++ )
            fprintf( fd, "%e ", kmatrix[ i*kN+j ] );
        fprintf( fd, "\n" );
    }
    fclose( fd );
}
*/

if ( this_task == 0 )
    fputs( "initialize kernel matrix ... done\n", stdout );
    if ( mean_flag )
        g2 = 1;
    else
        g2 = pow( g, 2 );
    for ( i=0; i<slice_num; i++ ) {
        for ( j=0; j<slice_index_num; j++ ) {
            if ( i == slice_index[j] ) {
                z1 = 0;
                while ( data[ z1*6+2 ] < i * dz ) z1++;
                z2 = z1;
                while ( data[ z2*6+2 ] < ( i+1 )*dz && z2<Particle[pt].num ) {
                    //fprintf( stdout, "%i: %f\n", z2, data[ z2*6+2 ] );
                    z2++;
                }
                z2--;
                if ( this_task == 0 ){
                    fprintf( stdout, "slice_info: z1=%li, z2=%li ( %.4lf ~ %.4lf )\n",
                            z1, z2, i*dz, (i+1)*dz );
                    for ( k=z1; k<z1+test_num; k++ ) {
                        fprintf( stdout, "( %15.5f %15.5f %15.5f ):   %e\n",
                                data[ k*6+0 ], data[ k*6+1 ], data[ k*6+2 ],
                                data[ k*6+3 ] );
                    }
                    for ( k=z2; k>z2-test_num; k-- ) {
                        fprintf( stdout, "( %15.5f %15.5f %15.5f ):   %e\n",
                                data[ k*6+0 ], data[ k*6+1 ], data[ k*6+2 ],
                                data[ k*6+3 ] );
                    }
                }
                for ( ii=0; ii<nxy[0]; ii++ )
                    for ( jj=0; jj<nxy[1]; jj++ ){
                        index = ii * nxy[1] + jj;
                        v1 [ index ] = 0;
                        }
                for ( ii=0; ii<anxy[0]; ii++ )
                    for ( jj=0; jj<anxy[1]; jj++ )
                        switch ( blk ) {
                            case IO_MAG:
                            case IO_MASS:
                            case IO_RHO:
                            case IO_VEL:
                                cgrid2.xg[ii][jj] = (ii * ag + slice_corner1[0]) / 1000;
                                cgrid2.yg[ii][jj] = (jj * ag + slice_corner1[1]) / 1000;
                                u[ii][jj] = 0;
                                v[ii][jj] = 0;
                        }
                for ( k=z1; k<z2; k++ ) {
                    if ( k % task_num == this_task ) {
                        index_i = ( long )( (data[ k*6+0 ]-slice_corner1[0]) / g );
                        index_j = ( long )( (data[ k*6+1 ]-slice_corner1[1]) / g );
                        index = index_i * nxy[1] + index_j;
                        if ( (proj_mode == 1 ) &&
                             ( (long)( ( data[ k*6+0 ]-slice_corner1[0]-h ) / g ) != index_i ||
                             (long)( ( data[ k*6+0 ]-slice_corner1[0]+h ) / g ) != index_i ||
                             (long)( ( data[ k*6+1 ]-slice_corner1[1]-h ) / g ) != index_j ||
                             (long)( ( data[ k*6+1 ]-slice_corner1[1]+h ) / g ) != index_j ) ) {
                            for ( ii=0; ii<kN; ii++ )
                                for ( jj=0; jj<kN; jj++ ) {
                                    index_i1 = (long) ( ( ( ii-kN/2 ) * h / kN + data[ k*6+0 ]-slice_corner1[0] ) / g );
                                    index_j1 = (long) ( ( ( jj-kN/2 ) * h / kN + data[ k*6+1 ]-slice_corner1[1] ) / g );
                                    if ( index_i1 < 0 || index_i1 > nxy[0] ||
                                         index_j1 < 0 || index_j1 > nxy[1] ) continue;
                                    index1 = index_i1 * nxy[1] + index_j1;
                                    v1[ index1 ] += data[ k*6+3 ] * kmatrix[ ii * kN +jj ] / g2;
                                }
                            //fprintf( stdout, "index1: %li\n", index1 );
                        }
                        else {
                            if ( index_i < 0 || index_i > nxy[0] ||
                                 index_j < 0 || index_j > nxy[1] ) continue;
                            v1[ index ] += data[ k*6+3 ] / g2;
#ifdef DEBUG
                            debug_l = index;
#endif
                        }
                        index_i = ( long )( (data[ k*6+0 ]-slice_corner1[0]) / ag );
                        index_j = ( long )( (data[ k*6+1 ]-slice_corner1[1]) / ag );
                        switch ( blk ) {
                            case IO_MAG:
                            case IO_MASS:
                            case IO_RHO:
                            case IO_VEL:
                                if ( index_i < 0 || index_i > anxy[0] ||
                                    index_j < 0 || index_j > anxy[1] ) continue;
                                u[index_i][index_j] += data[ k*6+4 ];
                                v[index_i][index_j] += data[ k*6+5 ];
                                break;
                        }
                    }
                }
                MPI_Reduce( v1, v2, nxy[0]*nxy[1], MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
if ( this_task == 0 ) {
               sprintf( fn_buf, "%s_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f",
                       out_picture_prefix,
                       slice_corner1[0], slice_corner2[0],
                       slice_corner1[1], slice_corner2[1],
                       i*dz, (i+1)*dz );
               fprintf( stdout, "picture name: %s\n", fn_buf );
               for ( ii=0; ii<nxy[0]; ii++ )
                   for ( jj=0; jj<nxy[1]; jj++ ){
                       v3[ii][jj]  = ( PLFLT ) ( v2[ ii*nxy[1]+jj ] );
                   }
               switch ( blk ) {
                    case IO_MAG:
                    case IO_RHO:
                    case IO_MASS:
                    case IO_VEL:
                        generate_2D_img( fn_buf, x_buf, "Mpc", v3, nxy, log_flag, xymm,
                            u, v, &cgrid2, anxy, flag_arrow );
                        break;
                    default:
                        generate_2D_img( fn_buf, x_buf, "Mpc", v3, nxy, log_flag, xymm,
                                NULL, NULL, NULL, NULL, 0 );
                        break;
               }
               if ( out_pic_data ) {
                    FILE *fd;
                    fd = fopen( fn_buf, "w" );
                    for ( ii=0; ii<nxy[0]; ii++ ) {
                        for ( jj=0; jj<nxy[1]; jj++ ) {
                           fprintf( fd, "%e ", v2[ ii*nxy[1]+jj ] );
                        }
                       fprintf( fd, "\n" );
                    }
                    fclose( fd );
               }
}
            }
        }
    }
    plFree2dGrid( v3, nxy[0], nxy[1] );
    free( v1 );
    if ( proj_mode == 1 )
        free( kmatrix );
if ( this_task == 0 )
        free( v2 );
if ( this_task == 0 )
    fputs( sep_str, stdout );
    free( data );
    switch ( blk ) {
        case IO_MASS:
        case IO_RHO:
            plFree2dGrid( u, anxy[0], anxy[1] );
            plFree2dGrid( v, anxy[0], anxy[1] );
            plFree2dGrid( cgrid2.xg, anxy[0], anxy[1] );
            plFree2dGrid( cgrid2.yg, anxy[0], anxy[1] );
            break;
        case IO_MAG:
        case IO_VEL:
            plFree2dGrid( u, anxy[0], anxy[1] );
            plFree2dGrid( v, anxy[0], anxy[1] );
            plFree2dGrid( cgrid2.xg, anxy[0], anxy[1] );
            plFree2dGrid( cgrid2.yg, anxy[0], anxy[1] );
        case IO_ACCEL:
            free( p );
            break;
    }
}

int compare_for_plot_position( const void *a, const void *b ){
    return ( ( (float*)b )[2] < ( ( float*)a )[2] ) ? 1 : 0;
}

void plot_position( int pt ) {
    PLFLT *x, *y, dz, xmin, xmax, ymin, ymax;
    PLINT i, j, k, num, array_len, index, s, test_num, z1, z2, z, flag;
    char fn_buf[50], buf[20];
    float *pos;
    fputs( sep_str, stdout );
    fprintf( stdout, "plot particle %i position ...\n", pt );
    pos = ( float* ) malloc( sizeof( float ) * Particle[pt].num * 3 );
    memcpy( pos, Particle[pt].pos, sizeof( float ) * Particle[pt].num *3 );
    test_num = 10;
    for ( i=0; i<test_num; i++ ) {
        fprintf( stdout, "%f %f %f\n", pos[i*3+0], pos[i*3+1], pos[i*3+2] );
    }
    if ( slice_num > 1 ){
        fputs( "Sorted by z: \n", stdout );
        qsort( (void*)pos, Particle[pt].num, sizeof( float )*3,
                &compare_for_plot_position );
        for ( i=0; i<test_num; i++ ) {
            fprintf( stdout, "%f %f %f\n", pos[i*3+0], pos[i*3+1], pos[i*3+2] );
        }
        fprintf( stdout, "\n" );
    }
    dz = header.BoxSize / slice_num;
    fprintf( stdout, "dz = %f\n", dz );
    for ( i=0; i<slice_num; i++ ) {
        flag = 0;
        if ( slice_index_num == 0 )
            flag = 1;
        else
        for ( j=0; j<=slice_index_num; j++ )
            if ( i == slice_index[j] ) {
                flag = 1;
                break;
            }
        if ( flag ) {
            fprintf( stdout, "plot slice: %i ( %.3f, %.3f )\n", i, i*dz, (i+1)*dz );
            z1 = 0;
            while ( pos[ z1*3+2 ] < i * dz ) z1++;
            z2 = z1;
            while ( pos[ z2*3+2 ] < ( i+1 )*dz && z2<Particle[pt].num ) z2++;
            z2--;
            fprintf( stdout, "slice_info: z1=%li, z2=%li\n", z1, z2 );
            for ( k=z1; k<z1+test_num; k++ ) {
                fprintf( stdout, "( %15.5f %15.5f %15.5f ):  \n",
                        pos[ k*3+0 ], pos[ k*3+1 ], pos[ k*3+2 ] );
            }
            for ( k=z2; k>z2-test_num; k-- ) {
                fprintf( stdout, "( %15.5f %15.5f %15.5f ):  \n",
                        pos[ k*3+0 ], pos[ k*3+1 ], pos[ k*3+2 ] );
            }
            sprintf( fn_buf, "%s_%i_%.2f_%i_%i.png", out_picture_prefix, pt, redshift, ( int )( i*dz ), ( int )((i+1)*dz) );
            plsdev( "pngcairo" );
            plsfnam( fn_buf );
            plinit();
            plcol0( 15 );
            xmin = 0;
            xmax = header.BoxSize;
            ymin = 0;
            ymax = header.BoxSize;
            plenv( xmin, xmax, ymin, ymax, 1, 0 );
            sprintf( buf, "z=%.2f", redshift );
            pllab( "", "", buf );
            num = z2 - z1;
            x = ( PLFLT* ) malloc( sizeof( PLFLT ) * num );
            y = ( PLFLT* ) malloc( sizeof( PLFLT ) * num );
            for ( k=z1; k<z2; k++ ){
                x[k-z1] = pos[ k*3+0 ];
                y[k-z1] = pos[ k*3+1 ];
                //fprintf( stdout, "( %f, %f )\n", x[k-z1], y[k-z1] );
            }
            plcol0( 2 );
            plssym( 0.35, 1 );
            plsym( num, x, y, 148 );
            plend();
            free( x );
            free( y );
        }
    }
    free( pos );
    fputs( sep_str, stdout );
}

void plot_box( unsigned int bit_flag, PLFLT *box ) {
    PLFLT x[2], y[2], z[2];
    if ( bit_flag & 1 ){
        x[0] = 0;
        y[0] = 0;
        z[0] = 0;
        x[1] = 0;
        y[1] = 0;
        z[1] = box[2];
        plline3( 2, x, y, z );
    }
    if ( bit_flag & 2 ){
        x[0] = 0;
        y[0] = box[1];
        z[0] = 0;
        x[1] = 0;
        y[1] = box[1];
        z[1] = box[2];
        plline3( 2, x, y, z );
    }
    if ( bit_flag & 4 ){
        x[0] = box[0];
        y[0] = 0;
        z[0] = 0;
        x[1] = box[0];
        y[1] = 0;
        z[1] = box[2];
        plline3( 2, x, y, z );
    }
    if ( bit_flag & 8 ){
        x[0] = box[0];
        y[0] = box[1];
        z[0] = 0;
        x[1] = box[0];
        y[1] = box[1];
        z[1] = box[2];
        plline3( 2, x, y, z );
    }
    if ( bit_flag & 16 ){
        x[0] = 0;
        y[0] = 0;
        z[0] = 0;
        x[1] = box[0];
        y[1] = 0;
        z[1] = 0;
        plline3( 2, x, y, z );
    }
    if ( bit_flag & 32 ){
        x[0] = 0;
        y[0] = 0;
        z[0] = box[2];
        x[1] = box[0];
        y[1] = 0;
        z[1] = box[2];
        plline3( 2, x, y, z );
    }
    if ( bit_flag & 64 ){
        x[0] = 0;
        y[0] = box[1];
        z[0] = 0;
        x[1] = box[0];
        y[1] = box[1];
        z[1] = 0;
        plline3( 2, x, y, z );
    }
    if ( bit_flag & 128 ){
        x[0] = 0;
        y[0] = box[1];
        z[0] = box[2];
        x[1] = box[0];
        y[1] = box[1];
        z[1] = box[2];
        plline3( 2, x, y, z );
    }
    if ( bit_flag & 256 ){
        x[0] = 0;
        y[0] = 0;
        z[0] = 0;
        x[1] = 0;
        y[1] = box[1];
        z[1] = 0;
        plline3( 2, x, y, z );
    }
    if ( bit_flag & 512 ){
        x[0] = 0;
        y[0] = 0;
        z[0] = box[2];
        x[1] = 0;
        y[1] = box[1];
        z[1] = box[2];
        plline3( 2, x, y, z );
    }
    if ( bit_flag & 1024 ){
        x[0] = box[0];
        y[0] = 0;
        z[0] = 0;
        x[1] = box[0];
        y[1] = box[1];
        z[1] = 0;
        plline3( 2, x, y, z );
    }
    if ( bit_flag & 2048 ){
        x[0] = box[0];
        y[0] = 0;
        z[0] = box[2];
        x[1] = box[0];
        y[1] = box[1];
        z[1] = box[2];
        plline3( 2, x, y, z );
    }
}

void plot_3d( long point_num, float *data, char *title_buf, char *x_buf, char *y_buf ) {
    char fn_buf[50], buf[20], title_buf2[300];
    unsigned int bit_flag, i;
    long num, n, png_index, naz, nal, iaz, ial;
    PLFLT box[3], *x, *y, *z, al_local, az_local;
if ( this_task == 0 )
    fputs( sep_str, stdout );
    box[0] = corner2[0] - corner1[0];
    box[1] = corner2[1] - corner1[1];
    box[2] = corner2[2] - corner1[2];
    num = 0;
    for ( i=0; i<point_num; i++ ) {
        if ( data[i*3+0] <= corner2[0] && data[i*3+0] >= corner1[0] &&
             data[i*3+1] <= corner2[1] && data[i*3+1] >= corner1[1] &&
             data[i*3+2] <= corner2[2] && data[i*3+2] >= corner1[2] )
            //fprintf( stdout, "%.2f %.2f %.2f\n", data[i*3+0], data[i*3+1], data[i*3+2] );
            num ++;
    }
    x = ( PLFLT* ) malloc( sizeof( PLFLT ) * num );
    y = ( PLFLT* ) malloc( sizeof( PLFLT ) * num );
    z = ( PLFLT* ) malloc( sizeof( PLFLT ) * num );
    if ( NULL == x || NULL == y || NULL == z ) {
if ( this_task == 0 )
        fprintf( stdout, "Failed to allocate x y z array!\n" );
        end_run( 7 );
    }
    n = 0;
    for ( i=0; i<point_num; i++ ) {
        if ( data[i*3+0] <= corner2[0] && data[i*3+0] >= corner1[0] &&
             data[i*3+1] <= corner2[1] && data[i*3+1] >= corner1[1] &&
             data[i*3+2] <= corner2[2] && data[i*3+2] >= corner1[2] ){
            x[n] = data[i*3+0] - corner1[0];
            y[n] = data[i*3+1] - corner1[1];
            z[n] = data[i*3+2] - corner1[2];
            n++;
        }
    }
if ( this_task == 0 )
    fprintf( stdout, "point number: %i\n", num );
    nal = ( long ) ( ( al[2]-al[0] ) / al[1] ) + 1;
    naz = ( long ) ( ( az[2]-az[0] ) / az[1] ) + 1;
if ( this_task == 0 )
    fprintf( stdout, "nal = %li, naz = %li\n", nal , naz );
sleep( 2 );
    for ( al_local=al[0]; al_local<=al[2]; al_local += al[1] ){
        for ( az_local=az[0]; az_local<=az[2]; az_local += az[1] ){
            ial = ( long ) ( ( al_local-al[0] ) / al[1] );
            iaz = ( long ) ( ( az_local-az[0] ) / al[1] );
            if ( ( ial * naz + iaz ) % task_num == this_task ) {
                plsdev( "pngcairo" );
                //sprintf( fn_buf, "%s_%i_%.2f_%.2f_%.2f.png", out_picture_prefix, pt, redshift, al_local, az_local );
                sprintf( fn_buf, "%s/%04i_%04i.png", out_picture_prefix, ial, iaz );
                fprintf( stdout, "process %3i plot: al_local=%8.2lf, az_local=%8.2lf, file name : %s\n", this_task,
                        al_local, az_local, fn_buf );
                plsfnam( fn_buf );
                plspage( 300, 300, 1024, 1024, 0, 0 );
                plinit();
                pladv( 0 );
                plcol0( 15 );
                plvpor( 0.0, 1.0, 0.0, 1.0 );
                //sprintf( title_buf2, "%s\n%s\n%s", title_buf, y_buf, x_buf );
                //pllab( x_buf, y_buf, title_buf );
                plwind( -box[0] * sqrt( 3.0 ),
                         box[0] * sqrt( 3.0 ),
                        -box[1] * sqrt( 3.0 ),
                         box[1] * sqrt( 3.0 ) );
                //plwind( -box[0], box[0], -box[1], box[1] );
                plw3d(  box[0], box[1], box[2],
                        0.0, box[0],
                        0.0, box[1],
                        0.0, box[2],
                        al_local, az_local );
                plbox3( "", "", 0.0, 0,
                        "", "", 0.0, 0,
                        "", "", 0.0, 0 );
                bit_flag = 1;
                for ( i=0; i<12; i++ ){
                    bit_flag = bit_flag << 1;
                    bit_flag ++;
                }
                plcol0( 15 );
                plwidth( 0.2 );
                plot_box( bit_flag, box );
                plcol0( 2 );
                plssym( 0.35, 1 );
                /*
                for ( i=0; i<num; i++ ) {
                    fprintf( stdout, "%f %f %f\n", x[i], y[i], z[i] );
                }
               */
                plpoin3( (PLINT)num, x, y, z, 1 );
                plend();
            }
        }
    }
    free( x );
    free( y );
    free( z );
if ( this_task == 0 )
    fputs( sep_str, stdout );
sleep( 2 );
}

void plot_3d_position( int pt ) {
    char title_buf[20], x_buf[100], y_buf[100];
    sprintf( title_buf, "z=%.2f", redshift );
    sprintf( x_buf, "boxsize(%.0f, %.0f, %.0f(Mpc))",
            (corner2[0]-corner1[0]) / 1000.0,
            (corner2[1]-corner1[1]) / 1000.0,
            (corner2[2]-corner1[2]) / 1000.0 );
    sprintf( y_buf, "" );
    plot_3d( Particle[pt].num, Particle[pt].pos, title_buf, x_buf, y_buf);
}

void plot_3d_scalar( int pt, enum iofields blk ) {
    float *density, *data, *p;
    char buf[20], title_buf[20], x_buf[100], y_buf[100];
    int i,j,k, tmp;
    long ii, num, n;
    srand( time(0) );
    density = ( float* ) malloc( sizeof( float ) * box[0] * box[1] * box[2] );
    memset( density, 0, sizeof(float) * box[0] * box[1] * box[2] );
    sprintf( title_buf, "z=%.2f", redshift );
    sprintf( x_buf, "boxsize(%.0f, %.0f, %.0f(Mpc))",
            (corner2[0]-corner1[0]) / 1000.0,
            (corner2[1]-corner1[1]) / 1000.0,
            (corner2[2]-corner1[2]) / 1000.0 );
    switch ( blk ) {
        case IO_U:
            get_dataset_name( blk, buf );
            if ( pt != 0 ) {
                get_dataset_name( blk, buf );
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            sprintf( y_buf, "%s(>%.2e)", buf, scalar_unit );
            p = Particle[pt].u;
            break;
        case IO_RHO:
            get_dataset_name( blk, buf );
            if ( pt != 0 ) {
                get_dataset_name( blk, buf );
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            sprintf( y_buf, "%s(>%.2e)", buf, scalar_unit );
            p = Particle[pt].rho;
            break;
        case IO_MASS:
            get_dataset_name( blk, buf );
            p = Particle[pt].m;
            sprintf( y_buf, "%s(>%.2e)", buf, scalar_unit );
            break;
        case IO_POT:
            get_dataset_name( blk, buf );
            p = Particle[pt].pot;
            sprintf( y_buf, "%s(>%.2e)", buf, scalar_unit );
            break;
        case IO_ELEC:
            get_dataset_name( blk, buf );
            if ( pt != 0 ) {
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            p = ( float* ) malloc( sizeof(float) * Particle[pt].num );
            for ( i=0; i<Particle[pt].num; i++ ) {
                p[i] = Particle[pt].rho[i] * Particle[pt].elec[i];
                //fprintf( stdout, "%e\n", p[i] );
            }
            sprintf( y_buf, "%s(>%.2e)", buf, scalar_unit );
            break;
        case IO_VEL:
            get_dataset_name( blk, buf );
            p = ( float* ) malloc( sizeof(float) * Particle[pt].num );
            for ( i=0; i<Particle[pt].num; i++ ) {
                p[i] = sqrt( pow( Particle[pt].vel[i*3+0], 2 ) +
                             pow( Particle[pt].vel[i*3+1], 2 ) +
                             pow( Particle[pt].vel[i*3+2], 2 ) );
            }
            sprintf( y_buf, "%s(>%.2e)", buf, scalar_unit );
            break;
        case IO_ACCEL:
            get_dataset_name( blk, buf );
            p = ( float* ) malloc( sizeof(float) * Particle[pt].num );
            for ( i=0; i<Particle[pt].num; i++ ) {
                p[i] = sqrt( pow( Particle[pt].accel[i*3+0], 2 ) +
                             pow( Particle[pt].accel[i*3+1], 2 ) +
                             pow( Particle[pt].accel[i*3+2], 2 ) );
            }
            sprintf( y_buf, "%s(>%.2e)", buf, scalar_unit );
            break;
        case IO_MAG:
            get_dataset_name( blk, buf );
            if ( pt != 0 ) {
                fprintf( stderr, "Particle %i hasn't field: \"%s\"\n", pt, buf );
                end_run( 3 );
            }
            p = ( float* ) malloc( sizeof(float) * Particle[pt].num );
            for ( i=0; i<Particle[pt].num; i++ ) {
                p[i] = sqrt( pow( Particle[pt].mag[i*3+0], 2 ) +
                             pow( Particle[pt].mag[i*3+1], 2 ) +
                             pow( Particle[pt].mag[i*3+2], 2 ) );
            }
            sprintf( y_buf, "%s(>%.2e)", buf, scalar_unit );
            break;
    }
    for ( ii=0; ii<Particle[pt].num; ii++ ) {
        i = ( long ) ( Particle[pt].pos[ii*3+0] / ( header.BoxSize / box[0] ) );
        j = ( long ) ( Particle[pt].pos[ii*3+1] / ( header.BoxSize / box[1] ) );
        k = ( long ) ( Particle[pt].pos[ii*3+2] / ( header.BoxSize / box[2] ) );
        density[ i*box[1]*box[2] + j*box[2] + k ] += ( p[ii] / scalar_unit );
    }
    num = 0;
    for ( i=0; i<box[0]; i++ )
        for ( j=0; j<box[1]; j++ )
            for ( k=0; k<box[2]; k++ ){
                num += (int)( density[ i*box[1]*box[2] + j*box[2] + k ] );
            }
if ( this_task == 0 )
    fprintf( stdout, "num = %li\n", num );
    data = ( float* ) malloc( sizeof( float ) * num * 3 );
    n = 0;
    for ( i=0; i<box[0]; i++ )
        for ( j=0; j<box[1]; j++ )
            for ( k=0; k<box[2]; k++ ) {
                tmp = ( int )(density[ i*box[1]*box[2] + j*box[2] + k ]);
                for ( ii=0; ii<tmp; ii++ ){
                    data[n*3+0] = rand()/(float)(RAND_MAX) * ( header.BoxSize / box[0] ) +
                                                            i * ( header.BoxSize / box[0] );
                    data[n*3+1] = rand()/(float)(RAND_MAX) * ( header.BoxSize / box[1] ) +
                                                            j * ( header.BoxSize / box[1] );
                    data[n*3+2] = rand()/(float)(RAND_MAX) * ( header.BoxSize / box[2] ) +
                                                            k * ( header.BoxSize / box[2] );
                    n++;
                }
            }
    plot_3d( num, data, title_buf, x_buf, y_buf );
    free( data );
    free( density );
    switch ( blk ) {
        case IO_MAG:
            if ( pt==0 ) free( p );
            break;
        case IO_ELEC:
            if ( pt==0 ) free( p );
            break;
        case IO_VEL:
        case IO_ACCEL:
            free( p );
            break;
    }
}

void plot_3d_multi( int flag ) {
    char fn_buf[50], buf[20], title_buf2[300];
    unsigned int bit_flag, i;
    long n, png_index, num[4], naz, nal, iaz, ial;
    PLFLT box[3], sx, sy, sz;
    PLFLT *x[4], *y[4], *z[4], az_local, al_local, *m;
if ( this_task == 0 )
    fputs( sep_str, stdout );
    box[0] = corner2[0] - corner1[0];
    box[1] = corner2[1] - corner1[1];
    box[2] = corner2[2] - corner1[2];
    switch ( flag ) {
        case 3:
            num[3] = 0;
            for ( i=0; i<TotNgroups; i++ ) {
                if ( group[i].Mass>1000 ) {
                    num[3] ++;
                }
            }
            x[3] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[3] );
            y[3] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[3] );
            z[3] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[3] );
            m = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[3] );
            n = 0;
            for ( i=0; i<TotNgroups; i++ ) {
                if ( group[i].Mass>1000 ) {
                    x[3][n] = group[i].CM[0];
                    y[3][n] = group[i].CM[1];
                    z[3][n] = group[i].CM[2];
                    m[n] = group[i].Mass;
                    n++;
                }
            }
if ( this_task == 0 )
            fprintf( stdout, "plot fof number: %i\n", num[3] );
        case 2:
            num[2] = 0;
            for ( i=0; i<Particle[1].num; i++ ) {
            if ( Particle[1].pos[i*3+0] <= corner2[0] && Particle[1].pos[i*3+0] >= corner1[0] &&
                Particle[1].pos[i*3+1] <= corner2[1] && Particle[1].pos[i*3+1] >= corner1[1] &&
                Particle[1].pos[i*3+2] <= corner2[2] && Particle[1].pos[i*3+2] >= corner1[2] )
                //fprintf( stdout, "%.2f %.2f %.2f\n", Particle[1].pos[i*3+0], Particle[1].pos[i*3+1], Particle[1].pos[i*3+2] );
                num[2] ++;
            }
            x[2] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[2] );
            y[2] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[2] );
            z[2] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[2] );
            n = 0;
            for ( i=0; i<Particle[1].num; i++ ) {
                if ( Particle[1].pos[i*3+0] <= corner2[0] && Particle[1].pos[i*3+0] >= corner1[0] &&
                    Particle[1].pos[i*3+1] <= corner2[1] && Particle[1].pos[i*3+1] >= corner1[1] &&
                    Particle[1].pos[i*3+2] <= corner2[2] && Particle[1].pos[i*3+2] >= corner1[2] ){
                    x[2][n] = Particle[1].pos[i*3+0] - corner1[0];
                    y[2][n] = Particle[1].pos[i*3+1] - corner1[1];
                    z[2][n] = Particle[1].pos[i*3+2] - corner1[2];
                    n++;
                }
            }
if ( this_task == 0 )
            fprintf( stdout, "plot halo number: %i\n", num[2] );
        case 1:
            num[1] = 0;
            for ( i=0; i<Particle[4].num; i++ ) {
            if ( Particle[4].pos[i*3+0] <= corner2[0] && Particle[4].pos[i*3+0] >= corner1[0] &&
                Particle[4].pos[i*3+1] <= corner2[1] && Particle[4].pos[i*3+1] >= corner1[1] &&
                Particle[4].pos[i*3+2] <= corner2[2] && Particle[4].pos[i*3+2] >= corner1[2] )
                //fprintf( stdout, "%.2f %.2f %.2f\n", Particle[4].pos[i*3+0], Particle[4].pos[i*3+1], Particle[4].pos[i*3+2] );
                num[1] ++;
            }
            x[1] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[1] );
            y[1] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[1] );
            z[1] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[1] );
            n = 0;
            for ( i=0; i<Particle[4].num; i++ ) {
                if ( Particle[4].pos[i*3+0] <= corner2[0] && Particle[4].pos[i*3+0] >= corner1[0] &&
                    Particle[4].pos[i*3+1] <= corner2[1] && Particle[4].pos[i*3+1] >= corner1[1] &&
                    Particle[4].pos[i*3+2] <= corner2[2] && Particle[4].pos[i*3+2] >= corner1[2] ){
                    x[1][n] = Particle[4].pos[i*3+0] - corner1[0];
                    y[1][n] = Particle[4].pos[i*3+1] - corner1[1];
                    z[1][n] = Particle[4].pos[i*3+2] - corner1[2];
                    n++;
                }
            }
if ( this_task == 0 )
            fprintf( stdout, "plot star number: %i\n", num[1] );
            num[0] = 0;
            for ( i=0; i<Particle[0].num; i++ ) {
                if ( Particle[0].pos[i*3+0] <= corner2[0] && Particle[0].pos[i*3+0] >= corner1[0] &&
                     Particle[0].pos[i*3+1] <= corner2[1] && Particle[0].pos[i*3+1] >= corner1[1] &&
                     Particle[0].pos[i*3+2] <= corner2[2] && Particle[0].pos[i*3+2] >= corner1[2] )
                    //fprintf( stdout, "%.2f %.2f %.2f\n", Particle[0].pos[i*3+0], Particle[0].pos[i*3+1], Particle[0].pos[i*3+2] );
                    num[0] ++;
            }
            x[0] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[0] );
            y[0] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[0] );
            z[0] = ( PLFLT* ) malloc( sizeof( PLFLT ) * num[0] );
            n = 0;
            for ( i=0; i<Particle[0].num; i++ ) {
            if ( Particle[0].pos[i*3+0] <= corner2[0] && Particle[0].pos[i*3+0] >= corner1[0] &&
                Particle[0].pos[i*3+1] <= corner2[1] && Particle[0].pos[i*3+1] >= corner1[1] &&
                Particle[0].pos[i*3+2] <= corner2[2] && Particle[0].pos[i*3+2] >= corner1[2] ){
                x[0][n] = Particle[0].pos[i*3+0] - corner1[0];
                y[0][n] = Particle[0].pos[i*3+1] - corner1[1];
                z[0][n] = Particle[0].pos[i*3+2] - corner1[2];
                n++;
                }
            }
if ( this_task == 0 )
            fprintf( stdout, "plot gas number: %i\n", num[0] );
            break;
    }
if ( this_task == 0 )
    fprintf( stdout, "point number: %i\n", num );
    nal = ( long ) ( ( al[2]-al[0] ) / al[1] ) + 1;
    naz = ( long ) ( ( az[2]-az[0] ) / az[1] ) + 1;
if ( this_task == 0 )
    fprintf( stdout, "nal = %li, naz = %li\n", nal , naz );
sleep( 2 );
    for ( al_local=al[0]; al_local<=al[2]; al_local += al[1] )
        for ( az_local=az[0]; az_local<=az[2]; az_local += az[1] ){
            ial = ( long ) ( ( al_local-al[0] ) / al[1] );
            iaz = ( long ) ( ( az_local-az[0] ) / al[1] );
            if ( ( ial * naz + iaz ) % task_num == this_task ) {
                plsdev( "pngcairo" );
                //sprintf( fn_buf, "%s_%i_%.2f_%.2f_%.2f.png", out_picture_prefix, pt, redshift, al_local, az_local );
                sprintf( fn_buf, "%s/%04i_%04i.png", out_picture_prefix, ial, iaz );
                fprintf( stdout, "process %3i plot: al_local=%8.2lf, az_local=%8.2lf, file name : %s\n", this_task,
                        al_local, az_local, fn_buf );
sleep( 2 );
                plsfnam( fn_buf );
                plspage( 300, 300, 1024, 1024, 0, 0 );
                plinit();
                pladv( 0 );
                plcol0( 15 );
                plvpor( 0.0, 1.0, 0.0, 1.0 );
                //sprintf( title_buf2, "%s\n%s\n%s", title_buf, y_buf, x_buf );
                //pllab( x_buf, y_buf, title_buf );
                plwind( -box[0] * sqrt( 3.0 ),
                         box[0] * sqrt( 3.0 ),
                        -box[1] * sqrt( 3.0 ),
                         box[1] * sqrt( 3.0 ) );
                //plwind( -box[0], box[0], -box[1], box[1] );
                /*
                plbox( "bcnt", 0.0, 0,
                        "bcnt", 0.0, 0 );
                        */
                plw3d(  box[0], box[1], box[2],
                        0.0, box[0],
                        0.0, box[1],
                        0.0, box[2],
                        al_local, az_local );
                plbox3( "", "", 0.0, 0,
                        "", "", 0.0, 0,
                        "", "", 0.0, 0 );
                bit_flag = 1;
                for ( i=0; i<12; i++ ){
                    bit_flag = bit_flag << 1;
                    bit_flag ++;
                }
                plcol0( 15 );
                plwidth( 0.2 );
                plot_box( bit_flag, box );
                switch ( flag ) {
                    case 3:
if ( this_task == 0 )
                        fputs( "plot fof ...\n", stdout );
                        plcol0( 15 );
                        plssym( 0.35, 20 );
                        plpoin3( (PLINT)(num[3]), x[3], y[3], z[3], 12 );
                        /*
                        for ( i=0; i<num[3]; i++ ){
                            sprintf( buf, "%.1f", m[i] );
                            sx = x[3][i];
                            sy = y[3][i];
                            sz = z[3][i];
                            plstring3( 1, &sx, &sy, &sz, buf );
                        }
                        */
                    case 2:
if ( this_task == 0 )
                        fputs( "plot halo ...\n", stdout );
                        plcol0( 9 );
                        plssym( 0.35, 1 );
                        plpoin3( (PLINT)(num[2]), x[2], y[2], z[2], 1 );
                    case 1:
if ( this_task == 0 )
                        fputs( "plot star ...\n", stdout );
                        plcol0( 2 );
                        plssym( 0.35, 1 );
                        plpoin3( (PLINT)(num[1]), x[1], y[1], z[1], 1 );
if ( this_task == 0 )
                        fputs( "plot gas ...\n", stdout );
                        plcol0( 1 );
                        plssym( 0.35, 1 );
                        plpoin3( (PLINT)(num[0]), x[0], y[0], z[0], 1 );
                }
                plend();
            }
        }
    switch ( flag ) {
        case 3:
            free( x[3] );
            free( y[3] );
            free( z[3] );
        case 2:
            free( x[2] );
            free( y[2] );
            free( z[2] );
        case 1:
            free( x[1] );
            free( y[1] );
            free( z[1] );
            free( x[0] );
            free( y[0] );
            free( z[0] );

    }
sleep( 2 );
if ( this_task == 0 )
    fputs( sep_str, stdout );
}

