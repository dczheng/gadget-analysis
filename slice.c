#include "allvars.h"

void slice_init() {

    int pt;
    long offset, num, index, i;
    ParticleData p_tmp;
    SphParticleData sphp_tmp;

    writelog( "determine slice info ...\n" );

    if ( End[proj_i] - Start[proj_i] !=
            End[proj_j] - Start[proj_j] )
        endruns( "Projection Region Must Be Square!..." );

    if ( End[0] == 0 ){
        Start[0] = 0;
        End[0] = BoxSize;
    }

    if ( End[1] == 0 ){
        Start[1] = 0;
        End[1] = BoxSize;
    }

    if ( End[2] == 0 ){
        Start[2] = 0;
        End[2] = BoxSize;
    }

    writelog( "StartX: %g, EndX: %g\n"
            "StartY: %g, EndY: %g\n"
            "StartZ: %g, EndZ: %g\n",
            Start[0], End[0],
            Start[1], End[1],
            Start[2], End[2] );

    for ( pt=0; pt<6; pt++ ) {

        offset = OffsetPart6[pt];
        num = NumPart6[pt];
        //printf( "%li %li\n", offset, num );
        index = offset;

        if ( num == 0 ) {
            SliceStart[pt] = -1;
            SliceEnd[pt] = -1;
            continue;
        }

        writelog( "particle %i: offset=%li, num=%li\n",
                pt, offset, num );

        for ( i=offset; i<offset+num; i++ ) {
            if ( P[i].Pos[0] >= Start[0] &&
                 P[i].Pos[0] <= End[0] &&
                 P[i].Pos[1] >= Start[1] &&
                 P[i].Pos[1] <= End[1] &&
                 P[i].Pos[2] >= Start[2] &&
                 P[i].Pos[2] <= End[2] ) {

                p_tmp = P[index];
                P[index] = P[i];
                P[i] = p_tmp;

                if ( pt == 0 ) {
                    sphp_tmp = SphP[ index ];
                    SphP[index] = SphP[i];
                    SphP[i] = sphp_tmp;
                }

                index ++;
            }
        }

        SliceStart[pt] = offset;
        SliceEnd[pt] = index;
    }
    writelog( "Slice Start: " );

    for ( pt=0; pt<6; pt++ ) {
        writelog( "%ld ", SliceStart[pt] );
    }

    writelog( "\n" );
    writelog( "Slice End: " );

    for ( pt=0; pt<6; pt++ ) {
        writelog( "%ld ", SliceEnd[pt] );
    }
    writelog( "\n" );

    writelog( "determine slice info ... done.\n" );
    put_sep0;

}

void make_slice_img( int pt, double *data, long NPart ) {

    double *img, *num, dx, dy, x, y, h, dh, lx, ly, v;
    int i, xi, yi, N, Nhalf, i1, i2, j1, j2, li, lj, PicSize;
    long start, end;

    PicSize = All.PicSize;
    writelog( "make slice imgage  ...\n" );

    reset_img();
    img = image.img;
    num = image.num;

    dx = dy = (End[proj_i] - Start[proj_i])/ PicSize;
    //writelog( "dx: %f, dy: %f\n", dx, dy  );

    N = All.KernelN;
    Nhalf = N / 2;
    h = SofteningTable[pt];
    dh = h / Nhalf;

    if ( NPart ) {
        start = 0;
        end = NPart;
    }
    else {
        start = SliceStart[pt];
        end = SliceEnd[pt];
    }

    for ( i=start; i<end; i++ ) {

        if ( NPart ) {
            x = data[i*3];
            y = data[i*3+1];
            v = data[i*3+2];
        }
        else {
            x = P[i].Pos[proj_i];
            y = P[i].Pos[proj_j];
            x -= Start[proj_i];
            y -= Start[proj_j];
            v = data[ i-SliceStart[pt] ];
        }

        //printf( "%g\n", v );
        xi = x / dx;
        yi = y / dy;
        //printf( "%i, %i\n", xi, yi );
        check_picture_index( xi );
        check_picture_index( yi );

        if ( All.KernelInterpolation == 0 ){
            img[ xi * PicSize + yi ] += v;
            num[ xi * PicSize + yi ] += 1;
            continue;
        }

        i1 = (int)(( x-h ) / dx);
        i2 = (int)(( x+h ) / dx);
        j1 = (int)(( y-h ) / dy);
        j2 = (int)(( y+h ) / dy);


        if ( i1 != xi || i2 != xi || j1 != yi || j2 != yi ) {
            /*
            writelog( "(%f, %f ), ( %i, %i), (%f, %f, %f, %f), (%i, %i, %i, %i)\n",
                x, y, xi, yi, x-h, x+h, y-h, y+h, i1, i2, j1, j2 );
                */
            for ( li=0; li<N; li++ )
                for ( lj=0; lj<N; lj++ ){
                        lx = x + ( li-Nhalf ) * dh;
                        ly = y + ( lj-Nhalf ) * dh;
                        i1 = lx / dx;
                        j1 = ly / dy;
                        if ( i1 < 0 || i1 >= PicSize ||
                                j1 < 0 || j1 >= PicSize ) continue;
                        img[ i1 * PicSize + j1 ] += v * KernelMat2D[pt][ li*N + lj ];
                        num[ i1 * PicSize + j1 ] += KernelMat2D[pt][ li*N + lj ];
                        //writelog( "%f %f\n", v, v * KernelMat2D[pt][li*N+lj] );
                }

        }
        else
            img[ xi * PicSize + yi ] += v;
            num[ xi * PicSize + yi ] += 1;
        //printf( "%g\n", img[ xi * PicSize + yi ] );

    }

    for ( i=0; i<SQR(All.PicSize); i++ )
        if ( num[i] != 0 )
            img[i] /= num[i];

    if ( All.UnitAreaSlice )
        for ( i=0; i<SQR(All.PicSize); i++ )
            img[i] /= SQR(dx);

    img_xmin = Start[ proj_i ];
    img_xmax = End[ proj_i ];
    img_ymin = Start[ proj_j ];
    img_ymax = End[ proj_j ];

}

void field_slice( int pt, double *data, char *name, long N ) {

    char buf[100];

    sprintf( buf, "`%s` slice\n", name );
    writelog( buf );

    sprintf( buf, "%s%s", OutputDir, name );
    create_dir( buf );
    sprintf( buf, "%s/%s_%03i.dat", buf, name, SnapIndex );

    if ( N )
        make_slice_img( 0, data, N );
    else
        make_slice_img( pt, data, 0 );

    write_img( buf, NULL );

}

void mag_slice() {
    int num, i;
    double *data;

    num = SliceEnd[0] - SliceStart[0];
    mymalloc2( data, sizeof( double ) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ ) {
        data[i] = get_B( i );
    }
    field_slice( 0, data, "MagneticField", 0 );
    myfree( data );
}

void mach_slice() {
    int num, i;
    double *data;

    num = SliceEnd[0] - SliceStart[0];
    mymalloc2( data, sizeof( double ) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ ) {
        data[i] = SphP[i].MachNumber;
    }
    field_slice( 0, data, "MachNumber", 0 );
    myfree( data );
}

void density_slice() {
    int num, i;
    long index;
    double *data, *data3;

    num = SliceEnd[0] - SliceStart[0];
    mymalloc2( data, sizeof( double ) * num );
    mymalloc2( data3, sizeof( double ) * num * 3 );

    for ( i=SliceStart[0]; i<num; i++ ) {
        data[i] = SphP[i].Density * ( g2c.g / CUBE( g2c.cm ) );
    }

    field_slice( 0, data, "Density", 0 );

    index = 0;
    for ( i=SliceStart[0]; i<num; i++ ) {
        if ( SphP[i].Temp >= 1e7 ) {
            data3[3*index] = P[i].Pos[proj_i] - Start[proj_i];
            data3[3*index+1] = P[i].Pos[proj_j] - Start[proj_j];
            data3[3*index+2] = SphP[i].Density * ( g2c.g / CUBE( g2c.cm ) );
            index ++;
        }
    }
    field_slice( 0, data3, "Density_hot", index );

    index = 0;
    for ( i=SliceStart[0]; i<num; i++ ) {
        if ( SphP[i].Temp < 1e7 && SphP[i].Temp >= 1e5 ) {
            data3[3*index] = P[i].Pos[proj_i] - Start[proj_i];
            data3[3*index+1] = P[i].Pos[proj_j] - Start[proj_j];
            data3[3*index+2] = SphP[i].Density * ( g2c.g / CUBE( g2c.cm ) );
            index ++;
        }
    }
    field_slice( 0, data3, "Density_warm-hot", index );

    index = 0;
    for ( i=SliceStart[0]; i<num; i++ ) {
        if ( ( SphP[i].Density / RhoBaryon ) < 1e3 && SphP[i].Temp < 1e5 ) {
            data3[3*index] = P[i].Pos[proj_i] - Start[proj_i];
            data3[3*index+1] = P[i].Pos[proj_j] - Start[proj_j];
            data3[3*index+2] = SphP[i].Density * ( g2c.g / CUBE( g2c.cm ) );
            index ++;
        }
    }
    field_slice( 0, data3, "Density_diffuse", index );

    index = 0;
    for ( i=SliceStart[0]; i<num; i++ ) {
        if ( ( SphP[i].Density / RhoBaryon ) >= 1e3 && SphP[i].Temp < 1e5 ) {
            data3[3*index] = P[i].Pos[proj_i] - Start[proj_i];
            data3[3*index+1] = P[i].Pos[proj_j] - Start[proj_j];
            data3[3*index+2] = SphP[i].Density * ( g2c.g / CUBE( g2c.cm ) );
            index ++;
        }
    }
    field_slice( 0, data3, "Density_condensed", index );

    myfree( data3 );
    myfree( data );
}

void temperature_slice() {
    int num, i;
    double *data;

    num = SliceEnd[0] - SliceStart[0];
    mymalloc2( data, sizeof( double ) * num );
    for ( i=SliceStart[0]; i<num; i++ ) {
        data[i] = SphP[i].Temp;
    }
    field_slice( 0, data, "Temperature", 0 );
    myfree( data );
}

void cren_slice() {
    int num, i;
    double *data;

    num = SliceEnd[0] - SliceStart[0];
    mymalloc2( data, sizeof( double ) * num );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ ) {
        data[i] = SphP[i].CRE_n * SphP[i].Density / guc.m_e / CUBE( g2c.cm );
    }
    field_slice( 0, data, "cre_n", 0 );
    myfree( data );
}

void radio_slice() {

    double dnu, *data, x, area, frac;
    int num, index, index1, i;
    char buf[20];

    dnu = log( All.NuMax/All.NuMin ) / ( All.NuNum-1 );
    num = SliceEnd[0] - SliceStart[0];

    x = log( All.Freq/All.NuMin ) / dnu;
    index = (int)(x);

    if ( index >= All.NuNum || index < 0 ) {
        endrun( 20190711 );
    }

    x -= index;

    area = (End[proj_i] - Start[proj_i]) / All.PicSize;
    area = SQR( area );

    index1 = ( index == All.NuNum-1 ) ? All.NuNum-1 : index + 1;

    mymalloc2( data, sizeof( double ) * num );

    frac = 1.0 / (4.0 * PI * SQR( LumDis * g2c.cm )) / ( area / SQR(ComDis) );
    for ( i=SliceStart[0]; i<SliceEnd[0]; i++ ) {
        data[i] = exp (
                log( PartRad[i*All.NuNum+index] ) * ( 1-x )
              + log( PartRad[i*All.NuNum+index1] ) * x
               ) * frac;
    }


    sprintf( buf, "radio_%.2f", All.Freq );
    field_slice( 0, data, buf, 0 );

    myfree( data );

}
