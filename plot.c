#include "allvars.h"

void generate_2D_img( char *fn_prefix, PLFLT **z, PLINT *nxy ){
    char buf[50];
    PLFLT zmin, zmax, colorbar_width, colorbar_heigh;
    PLINT i,j;
    sprintf( buf, "%s.png", fn_prefix );
    plsdev( "png" );
    plsfnam( buf );
    sprintf( buf, "z=%.2f", redshift );
    plinit();
    plcol0( 2 );
    plenv( 0.0, (PLFLT)header.BoxSize/1000 ,
           0.0, (PLFLT)header.BoxSize/1000 ,
           1, 0 );
    pllab( "Mpc", "", buf );
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
    for ( i=0; i<nxy[0]; i++ )
        for ( j=0; j<nxy[1]; j++ ){
            if ( z[i][j] == 0 ) z[i][j] = zmin;
            z[i][j] = log10( z[i][j] );
        }
    zmin = log10( zmin );
    zmax = log10( zmax );
    fprintf( stdout, "zmin=%e, zmax=%e\n", zmin, zmax );
    plimage( (PLFLT_MATRIX)z, nxy[0], nxy[1],
            0.0, ( PLFLT ) header.BoxSize/1000,
            0.0, ( PLFLT ) header.BoxSize/1000,
            zmin, zmax,
            0.0, ( PLFLT ) header.BoxSize/1000,
            0.0, ( PLFLT ) header.BoxSize/1000 );
    //set colorbar
    PLINT n_axis = 1;
    PLCHAR_VECTOR axis_opts[] = { "bcvtm", };
    PLFLT axis_ticks[1] = { 0.0, };
    PLINT axis_subticks[1] = { 0, };
    PLINT num_values[1] = { 10 };
    PLFLT *values[1];
    values[0] = ( PLFLT* ) malloc( num_values[0] * sizeof( PLFLT ) );
    for ( i=0; i<num_values[0]; i++ ) {
        values[0][i] = i * ( zmax-zmin ) / num_values[0] + zmin;
    }
    plcolorbar( &colorbar_width, &colorbar_heigh,
            PL_COLORBAR_IMAGE, 0,
            0.05, 0, 0.05, 0.9, 0, 1, 1, 0.0, 0.0, 0.0, 0.0,
            0, NULL, NULL,
            n_axis, axis_opts,
            axis_ticks, axis_subticks,
            num_values, (PLFLT_MATRIX) values );
    free( values[0] );
    //set colorbar
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

double z_integrand( double z, void *params ) {
    double r, *p;
    p = ( double* ) params;
    r = sqrt( pow( p[7] - p[2], 2 ) +
              pow( p[8] - p[3], 2 ) +
              pow( z, 2 ) );
              //pow( p[4], 2 ) );
    return kernel( z, p[6] );
}

double y_integrand( double y, void *params ) {
    double epsabs, epsrel, abserr, result, *p;
    size_t subinter;
    epsabs = epsrel = 1e-7;
    subinter = 10000;
    p = ( double* ) params;
    gsl_function F;
    gsl_integration_workspace *integration_workspace =
        gsl_integration_workspace_alloc( subinter );
    p[8] = y;
    F.function = &z_integrand;
    F.params = params;
    gsl_integration_qag( &F, -p[6], p[6], epsabs, epsrel,
            subinter, GSL_INTEG_GAUSS31, integration_workspace, &result, &abserr );
    /*
    gsl_integration_qag( &F, 0, p[6], epsabs, epsrel,
            subinter, GSL_INTEG_GAUSS31, integration_workspace, &result, &abserr );
    fprintf( stdout, "result: %lf\n", result );
            */
    gsl_integration_workspace_free( integration_workspace );
    return result;
}

double x_integrand( double x, void *params ) {
    double epsabs, epsrel, abserr, result, *p;
    size_t subinter;
    epsabs = epsrel = 1e-7;
    subinter = 10000;
    gsl_function F;
    p = ( double* ) params;
    p[7] = x;
    F.function = &y_integrand;
    F.params = params;
    gsl_integration_workspace *integration_workspace =
        gsl_integration_workspace_alloc( subinter );
    gsl_integration_qag( &F, p[1]-p[5]/2.0, p[1]+p[5]/2.0, epsabs, epsrel,
            subinter, GSL_INTEG_GAUSS31, integration_workspace, &result, &abserr );
    gsl_integration_workspace_free( integration_workspace );
    //fprintf( stdout, "y: %e\n", abserr );
    return result;
}

double los_integration( double *params ) {
    double epsabs, epsrel, abserr, result;
    size_t subinter;
    epsabs = epsrel = 1e-7;
    subinter = 10000;
    gsl_function F;
    F.function = &x_integrand;
    F.params = ( void* )params;
    gsl_integration_workspace *integration_workspace =
        gsl_integration_workspace_alloc( subinter );
    gsl_integration_qag( &F, params[0]-params[5]/2.0, params[0]+params[5]/2.0, epsabs, epsrel,
            subinter, GSL_INTEG_GAUSS31, integration_workspace, &result, &abserr );
    gsl_integration_workspace_free( integration_workspace );
    //fprintf( stdout, "x: %e\n", abserr );
    return result;
}

void plot_slice( int pt, enum iofields blk ){
    float *data, dz, r, *p;
    char fn_buf[50], buf[20];
    long i, N, j, z1, z2, k, test_num, ii, jj, index, index2;
    double g, h, params[9], rr1, rr2, rr3, rr4;
    PLINT nxy[2];
    PLFLT **v3, xymm[4];
    double *v1, *v2;;
    /*
    struct ngb_struct {
        int N, *ngb, len;
    } *ngb;
    */

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
            p = Particle[pt].pot;
            break;
        case IO_MN:
            p = Particle[pt].mn;
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
    plAlloc2dGrid( &v3, nxy[0], nxy[1] );
    v1 = ( double * ) malloc( sizeof( double ) * nxy[0] * nxy[1] );
    if ( this_task == 0 )
        v2 = ( double * ) malloc( sizeof( double ) * nxy[0] * nxy[1] );
    N = Particle[0].num;
    if ( this_task == 0 )
    fputs( sep_str, stdout );
    get_dataset_name( blk, buf );
    if ( this_task == 0 )
    fprintf( stdout, "plot partilcle %i field: \"%s\" ... \n", pt, buf );
    data = ( float* ) malloc( sizeof( float ) * N * 4 );
    for ( i=0; i<N; i++ ){
        data[ i*4+0 ] = Particle[pt].pos[ i*3+0 ];
        data[ i*4+1 ] = Particle[pt].pos[ i*3+1 ];
        data[ i*4+2 ] = Particle[pt].pos[ i*3+2 ];
        data[ i*4+3 ] = p[ i ];
    }
    test_num = 10;
    if ( this_task == 0 ){
        fputs( sep_str, stdout );
        for ( i=0; i<test_num; i++ ) {
            fprintf( stdout, "( %15.5f %15.5f %15.5f ):   %e\n",
                    data[ i*4+0 ], data[ i*4+1 ], data[ i*4+2 ],
                    data[ i*4+3 ] );
        }
    }
    qsort( (void*)data, N, sizeof( float )*4,
                &compare_for_plot_slice );
    if ( this_task == 0 ){
        fputs( "Sorted by z: \n", stdout );
        for ( i=0; i<test_num; i++ ) {
            fprintf( stdout, "( %15.5f %15.5f %15.5f ):   %e\n",
                    data[ i*4+0 ], data[ i*4+1 ], data[ i*4+2 ],
                    data[ i*4+3 ] );
        }
    }
    if ( this_task == 0 )
        fputs( sep_str, stdout );
    dz = header.BoxSize / slice_num;
    if ( this_task == 0 )
        fprintf( stdout, "slice info: dz=%f\n", dz );
    g = header.BoxSize / nxy[0];
    h = 2.5;
    switch( blk ) {
        case IO_RHO:
            for ( i=0; i<N; i++ )
                data[ i*4+3 ] *= g * g * h;
            break;
    }
    if ( this_task == 0 )
        fprintf( stdout, "g=%f, h=%f\n", g, h );
    params[5] = g;
    params[6] = h;
    /*
    ngb = ( struct ngb_struct* ) malloc( sizeof ( struct ngb_struct ) * pic_xsize * pic_ysize );
    for ( i=0; i<pic_xsize*pic_ysize; i++ ){
        ngb[i].N = 0;
        ngb[i].len = 0;
    }
    */
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
                if ( this_task == 0 ){
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
                }
                /*
                for ( ii=0; ii< nxy[0]; ii++ )
                    for ( jj=0; jj<nxy[1]; jj++ )
                        v3[ii][jj] = 0;
                for ( k=z2; k>z1; k-- ) {
                    ii = ( PLINT )(data[ k*4+0 ] / ( header.BoxSize / nxy[0] ));
                    jj = ( PLINT )(data[ k*4+1 ] / ( header.BoxSize / nxy[1] ));
                    v3[ii][jj] += ( PLFLT ) ( data[ k*4+3 ]  );
                    //fprintf( stdout, "(i=%i, j=%i): %f\n", ii, jj, v3[ii][jj] );
                    //r = ( dz - (data[ k*4+2 ] - i*dz) ) / dz;
                    //r = ( float )( z2-k ) / ( z2-z1-1 ) ;
                    //v3[ii][jj] += ( PLFLT ) ( data[ k*4+3 ] * pow( r,1 ) );
                    //fprintf( stdout, "%f\n", r );
                }
                */
                /*
                for ( ii=0; ii<nxy[0]; ii++ )
                    for ( jj=0; jj<nxy[1]; jj++ )
                        v1[ii * nxy[1] + jj] = 0;
                for ( ii=0; ii<nxy[0]; ii++ )
                    for ( jj=0; jj<nxy[1]; jj++ ) {
                        index = ii*nxy[1] + jj;
                        for ( k=z1; k<z2; k++ ) {
                            if ( ngb[ index ].len == 0 ) {
                                ngb[ index ].len = 32;
                                ngb[ index ].ngb = ( int* ) malloc( sizeof( int ) *
                                        ngb[ index ].len );
                            }
                            if ( ngb[ index ].N == ngb[ index ].len ) {
                                ngb[ index ].len += 32;
                                ngb[ index ].ngb = ( int* ) realloc( ngb[ index ].ngb,
                                        sizeof( int ) * ngb[ index ].len);
                            }
                            params[0] = ii * g;
                            params[1] = jj * g;
                            rr1 = sqrt( pow( data[k*4+0]-params[0], 2 ) +
                                      pow( data[k*4+1]-params[1], 2 ) +
                                      pow( data[k*4+2], 2 ));
                            rr2 = sqrt( pow( data[k*4+0]-params[0]-g, 2 ) +
                                      pow( data[k*4+1]-params[1], 2 ) +
                                      pow( data[k*4+2], 2 ));
                            rr3 = sqrt( pow( data[k*4+0]-params[0], 2 ) +
                                      pow( data[k*4+1]-params[1]-g, 2 ) +
                                      pow( data[k*4+2], 2 ));
                            rr4 = sqrt( pow( data[k*4+0]-params[0]-g, 2 ) +
                                      pow( data[k*4+1]-params[1]-g, 2 ) +
                                      pow( data[k*4+2], 2 ));
                            if  ( ( ( data[k*4+0] < ii * g - g ) ||
                                 ( data[k*4+0] > ii * g + g ) ||
                                 ( data[k*4+1] < jj * g - g ) ||
                                 ( data[k*4+1] > jj * g + g ) ) &&
                                    kernel( rr1, h ) == 0 &&
                                    kernel( rr2, h ) == 0 &&
                                    kernel( rr3, h ) == 0 &&
                                    kernel( rr4, h ) == 0 )
                                 continue;
                            ngb[ index ].ngb[ ngb[ index ].N ] = k;
                            ngb[ index ].N++;
                        }
                 //       fprintf( stdout, "%i %i %i\n", ii, jj, ngb[index].N );
                    }
                for ( ii=0; ii<nxy[0]; ii++ ){
                    for ( jj=0; jj<nxy[1]; jj++ ) {
                        index = ii * nxy[1] + jj;
                        if ( index % task_num == this_task ){
                            params[0] = ii * g;
                            params[1] = jj * g;
                            for ( k=0; k<ngb[index].N; k++ ) {
                                index2 = ngb[index].ngb[k];
                                params[2] = data[ index2*4+0 ];
                                params[3] = data[ index2*4+1 ];
                                params[4] = data[ index2*4+2 ];
                                v1[ii * nxy[1] + jj] +=  pow( g, -2 ) * pow( h, -3 ) *
                                    data[ k*4+3 ] * los_integration( (void*) params );
                                //fprintf( stdout, "%i %i %e %e %e\n", ii, jj, tmp, v1[ ii * nxy[1] + jj ], data[ k*4+3 ] );
                            }
                        }
                    }
                }
                for ( ii=0; ii<nxy[0]; ii++ )
                    for ( jj=0; jj<nxy[1]; jj++ ) {
                        index = ii*nxy[1] + jj;
                        free( ngb[index].ngb );
                    }
                    */
                for ( ii=0; ii<nxy[0]; ii++ ){
                    for ( jj=0; jj<nxy[1]; jj++ ) {
                        index = ii * nxy[1] + jj;
                        if ( index % task_num == this_task ){
                            params[0] = ii * g;
                            params[1] = jj * g;
                            for ( k=z1; k<z2; k++ ) {
                            rr1 = sqrt( pow( data[k*4+0]-params[0], 2 ) +
                                      pow( data[k*4+1]-params[1], 2 ) +
                                      pow( data[k*4+2], 2 ));
                            rr2 = sqrt( pow( data[k*4+0]-params[0]-g, 2 ) +
                                      pow( data[k*4+1]-params[1], 2 ) +
                                      pow( data[k*4+2], 2 ));
                            rr3 = sqrt( pow( data[k*4+0]-params[0], 2 ) +
                                      pow( data[k*4+1]-params[1]-g, 2 ) +
                                      pow( data[k*4+2], 2 ));
                            rr4 = sqrt( pow( data[k*4+0]-params[0]-g, 2 ) +
                                      pow( data[k*4+1]-params[1]-g, 2 ) +
                                      pow( data[k*4+2], 2 ));
                            if  ( ( ( data[k*4+0] < ii * g - g ) ||
                                 ( data[k*4+0] > ii * g + g ) ||
                                 ( data[k*4+1] < jj * g - g ) ||
                                 ( data[k*4+1] > jj * g + g ) ) &&
                                    kernel( rr1, h ) == 0 &&
                                    kernel( rr2, h ) == 0 &&
                                    kernel( rr3, h ) == 0 &&
                                    kernel( rr4, h ) == 0 )
                                 continue;
                                params[2] = data[ k*4+0 ];
                                params[3] = data[ k*4+1 ];
                                params[4] = data[ k*4+2 ];
                                v1[ index ] +=  pow( g, -2 ) *
                                    data[ k*4+3 ] * los_integration( (void*) params );
                                //fprintf( stdout, "%i %i %e %e %e\n", ii, jj, tmp, v1[ ii * nxy[1] + jj ], data[ k*4+3 ] );
                            }
                        }
                    }
                }
                MPI_Reduce( v1, v2, nxy[0]*nxy[1], MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
                if ( this_task == 0 ) {
                    sprintf( fn_buf, "%s%i_%i", out_picture_prefix, ( int )( i*dz ), ( int )((i+1)*dz) );
                    for ( ii=0; ii<nxy[0]; ii++ )
                        for ( jj=0; jj<nxy[1]; jj++ ){
                            v3[ii][jj]  = ( PLFLT ) ( v2[ ii*nxy[1]+jj ] );
                        }
                    generate_2D_img( fn_buf, v3, nxy );
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
    plFree2dGrid( v3, nxy[0], nxy[1] );
    //free( ngb );
    free( v1 );
    if ( this_task == 0 )
        free( v2 );
    if ( this_task == 0 )
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
    fputs( "Sorted by z: \n", stdout );
    qsort( (void*)pos, Particle[pt].num, sizeof( float )*3,
                &compare_for_plot_position );
    for ( i=0; i<test_num; i++ ) {
        fprintf( stdout, "%f %f %f\n", pos[i*3+0], pos[i*3+1], pos[i*3+2] );
    }
    fprintf( stdout, "\n" );
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
    long num, n, png_index;
    PLFLT box[3], *x, *y, *z, al_local, az_local;
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
    fprintf( stdout, "point number: %i\n", num );
    png_index = 0;
    for ( al_local=al[0]; al_local<=al[2]; al_local += al[1] ){
        for ( az_local=az[0]; az_local<=az[2]; az_local += az[1] ){
            plsdev( "pngcairo" );
            //sprintf( fn_buf, "%s_%i_%.2f_%.2f_%.2f.png", out_picture_prefix, pt, redshift, al_local, az_local );
            sprintf( fn_buf, "%s/%04i.png", out_picture_prefix, png_index );
            png_index++;
            fprintf( stdout, "plot: al_local=%lf, az_local=%lf, file name : %s\n",
                    al_local, az_local, fn_buf );
            plsfnam( fn_buf );
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
    free( x );
    free( y );
    free( z );
    fputs( sep_str, stdout );
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
    long n, png_index, num[4];
    PLFLT box[3], sx, sy, sz;
    PLFLT *x[4], *y[4], *z[4], az_local, al_local, *m;
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
            fprintf( stdout, "plot gas number: %i\n", num[0] );
            break;
    }
    png_index = 0;
    for ( al_local=al[0]; al_local<=al[2]; al_local += al[1] )
        for ( az_local=az[0]; az_local<=az[2]; az_local += az[1] ){
            plsdev( "pngcairo" );
            //sprintf( fn_buf, "%s_%i_%.2f_%.2f_%.2f.png", out_picture_prefix, pt, redshift, al_local, az_local );
            sprintf( fn_buf, "%s/%04i.png", out_picture_prefix, png_index );
            png_index++;
            fprintf( stdout, "plot: al_local=%lf, az_local=%lf, file name : %s\n",
                    al_local, az_local, fn_buf );
            plsfnam( fn_buf );
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
                    fputs( "plot halo ...\n", stdout );
                    plcol0( 9 );
                    plssym( 0.35, 1 );
                    plpoin3( (PLINT)(num[2]), x[2], y[2], z[2], 1 );
                case 1:
                    fputs( "plot star ...\n", stdout );
                    plcol0( 2 );
                    plssym( 0.35, 1 );
                    plpoin3( (PLINT)(num[1]), x[1], y[1], z[1], 1 );
                    fputs( "plot gas ...\n", stdout );
                    plcol0( 1 );
                    plssym( 0.35, 1 );
                    plpoin3( (PLINT)(num[0]), x[0], y[0], z[0], 1 );
            }
            plend();
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
    fputs( sep_str, stdout );
}

