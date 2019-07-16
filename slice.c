#include "allvars.h"

void slice_init() {

    int pt;
    long offset, num, index, i;
    struct Particle_Data p_tmp;
    struct Sph_Particle_Data sphp_tmp;

    writelog( "determine slice info ...\n" );

    if ( All.End[All.proj_i] - All.Start[All.proj_i] !=
            All.End[All.proj_j] - All.Start[All.proj_j] )
        endruns( "Projection Region Must Be Square!..." );

    if ( All.End[0] == 0 ){
        All.Start[0] = 0;
        All.End[0] = All.BoxSize;
    }

    if ( All.End[1] == 0 ){
        All.Start[1] = 0;
        All.End[1] = All.BoxSize;
    }

    if ( All.End[2] == 0 ){
        All.Start[2] = 0;
        All.End[2] = All.BoxSize;
    }

    writelog( "StartX: %g, EndX: %g\n"
            "StartY: %g, EndY: %g\n"
            "StartZ: %g, EndZ: %g\n",
            All.Start[0], All.End[0],
            All.Start[1], All.End[1],
            All.Start[2], All.End[2] );

    for ( pt=0; pt<6; pt++ ) {

        offset = OffsetPart6[pt];
        num = NumPart6[pt];
        //printf( "%li %li\n", offset, num );
        index = offset;

        if ( num == 0 ) {
            All.SliceStart[pt] = -1;
            All.SliceEnd[pt] = -1;
            continue;
        }

        writelog( "particle %i: offset=%li, num=%li\n",
                pt, offset, num );

        for ( i=offset; i<offset+num; i++ ) {
            if ( P[i].Pos[0] >= All.Start[0] &&
                 P[i].Pos[0] <= All.End[0] &&
                 P[i].Pos[1] >= All.Start[1] &&
                 P[i].Pos[1] <= All.End[1] &&
                 P[i].Pos[2] >= All.Start[2] &&
                 P[i].Pos[2] <= All.End[2] ) {

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

        All.SliceStart[pt] = offset;
        All.SliceEnd[pt] = index;
    }
    writelog( "Slice Start: " );

    for ( pt=0; pt<6; pt++ ) {
        writelog( "%ld ", All.SliceStart[pt] );
    }

    writelog( "\n" );
    writelog( "Slice End: " );

    for ( pt=0; pt<6; pt++ ) {
        writelog( "%ld ", All.SliceEnd[pt] );
    }
    writelog( "\n" );

    writelog( "determine slice info ... done.\n" );
    put_sep0;

}

void make_slice_img( int pt, double *data ) {

    double *img, *num, dx, dy, x, y, h, dh, lx, ly, v;
    int i, xi, yi, N, Nhalf, i1, i2, j1, j2, li, lj, PicSize;

    PicSize = All.PicSize;
    writelog( "make slice imgage  ...\n" );

    reset_img();
    img = image.img;
    num = image.num;

    dx = dy = (All.End[All.proj_i] - All.Start[All.proj_i])/ PicSize;
    //writelog( "dx: %f, dy: %f\n", dx, dy  );

    N = All.KernelN;
    Nhalf = N / 2;
    h = All.SofteningTable[pt];
    dh = h / Nhalf;

    for ( i=All.SliceStart[pt]; i<All.SliceEnd[pt]; i++ ) {

        x = P[i].Pos[All.proj_i];
        y = P[i].Pos[All.proj_j];
        x -= All.Start[All.proj_i];
        y -= All.Start[All.proj_j];
        v = data[ i-All.SliceStart[pt] ];

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
                        img[ i1 * PicSize + j1 ] += v * All.KernelMat2D[pt][ li*N + lj ];
                        num[ i1 * PicSize + j1 ] += All.KernelMat2D[pt][ li*N + lj ];
                        //writelog( "%f %f\n", v, v * All.KernelMat2D[pt][li*N+lj] );
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

    img_xmin = All.Start[ All.proj_i ];
    img_xmax = All.End[ All.proj_i ];
    img_ymin = All.Start[ All.proj_j ];
    img_ymax = All.End[ All.proj_j ];

}

void field_slice( int pt, double *data, char *name ) {

    char buf[100];

    sprintf( buf, "`%s` slice\n", name );
    writelog( buf );

    sprintf( buf, "%s%s", All.OutputDir, name );
    create_dir( buf );
    sprintf( buf, "%s/%s_%03i.dat", buf, name, All.SnapIndex );

    make_slice_img( pt, data );

    write_img1( buf, NULL );

}

void mag_slice() {
    int num, i;
    double *data;

    num = All.SliceEnd[0] - All.SliceStart[0];
    mymalloc2( data, sizeof( double ) * num );
    for ( i=All.SliceStart[0]; i<All.SliceEnd[0]; i++ ) {
        data[i] = get_B( i );
    }
    field_slice( 0, data, "MagneticField" );
    myfree( data );
}

void mach_slice() {
    int num, i;
    double *data;

    num = All.SliceEnd[0] - All.SliceStart[0];
    mymalloc2( data, sizeof( double ) * num );
    for ( i=All.SliceStart[0]; i<All.SliceEnd[0]; i++ ) {
        data[i] = SphP[i].MachNumber;
    }
    field_slice( 0, data, "MachNumber" );
    myfree( data );
}

void density_slice() {
    int num, i;
    double *data;

    num = All.SliceEnd[0] - All.SliceStart[0];
    mymalloc2( data, sizeof( double ) * num );
    for ( i=All.SliceStart[0]; i<num; i++ ) {
        data[i] = SphP[i].Density * ( g2c.g / CUBE( g2c.cm ) );
    }

    field_slice( 0, data, "Density" );
    myfree( data );
}

void temperature_slice() {
    int num, i;
    double *data;

    num = All.SliceEnd[0] - All.SliceStart[0];
    mymalloc2( data, sizeof( double ) * num );
    for ( i=All.SliceStart[0]; i<num; i++ ) {
        data[i] = SphP[i].Temp;
    }
    field_slice( 0, data, "Temperature" );
    myfree( data );
}

void cren_slice() {
    int num, i;
    double *data;

    num = All.SliceEnd[0] - All.SliceStart[0];
    mymalloc2( data, sizeof( double ) * num );
    for ( i=All.SliceStart[0]; i<All.SliceEnd[0]; i++ ) {
        data[i] = SphP[i].CRE_n * SphP[i].Density / guc.m_e / CUBE( g2c.cm );
    }
    field_slice( 0, data, "cre_n" );
    myfree( data );
}

void radio_slice() {

    double dnu, *data, x, area, frac;
    int num, index, index1, i;
    char buf[20];

    dnu = log( All.NuMax/All.NuMin ) / ( All.NuNum-1 );
    num = All.SliceEnd[0] - All.SliceStart[0];

    x = log( All.Freq/All.NuMin ) / dnu;
    index = (int)(x);

    if ( index >= All.NuNum || index < 0 ) {
        endrun( 20190711 );
    }

    x -= index;

    area = (All.End[All.proj_i] - All.Start[All.proj_i]) / All.PicSize;
    area = SQR( area );

    index1 = ( index == All.NuNum-1 ) ? All.NuNum-1 : index + 1;

    mymalloc2( data, sizeof( double ) * num );

    frac = 1.0 / (4.0 * PI * SQR( All.LumDis * g2c.cm )) / ( area / SQR(All.ComDis) );
    for ( i=All.SliceStart[0]; i<All.SliceEnd[0]; i++ ) {
        data[i] = exp (
                log( PartRad[i*All.NuNum+index] ) * ( 1-x )
              + log( PartRad[i*All.NuNum+index1] ) * x
               ) * frac;
    }


    sprintf( buf, "radio_%.2f", All.Freq );
    field_slice( 0, data, buf );

    myfree( data );

}
