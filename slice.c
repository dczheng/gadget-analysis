#include "allvars.h"

void slice() {
    char buf[200];
    int pt;
    long long offset, num, index, i;
    struct particle_data p_tmp;
    struct sph_particle_data sphp_tmp;
    writelog( "determine slice info ...\n" );

    if ( All.End[proj_i] - All.Start[proj_i] !=
            All.End[proj_j] - All.Start[proj_j] ) {
        printf( "Projection Region Must Be Square!...\n" );
        endrun( 20180424 );
    }

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
        offset = get_particle_offset( pt );
        num = get_particle_num( pt );
        index = offset;
        if ( num == 0 ) {
            SliceStart[pt] = -1;
            SliceEnd[pt] = -1;
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
                P[i] = P[index];
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
}

void make_slice_img( int pt ) {
    double *img, dx, dy, x, y, h, dh, lx, ly, v;
    int i, j, xi, yi, N, Nhalf, i1, i2, j1, j2, li, lj, PicSize;
    PicSize = All.PicSize;
    writelog( "make slice imgage  ...\n" );
    img = image.img;
    memset( img, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = (All.End[proj_i] - All.Start[proj_i])/ PicSize;
    image.DataMin = image.ImgMin = DBL_MAX;
    image.DataMax = image.ImgMax = DBL_MIN;
    N = All.KernelN;
    Nhalf = N / 2;
    h = All.SofteningTable[pt];
    dh = h / Nhalf;
    for ( i=SliceStart[pt]; i<SliceEnd[pt]; i++ ){
        x = P[i].Pos[proj_i];
        y = P[i].Pos[proj_j];
        x -= All.Start[proj_i];
        y -= All.Start[proj_j];
        v = image.data[ i-SliceStart[pt] ];
        //printf( "%g\n", v );
        xi = x / dx;
        yi = y / dy;
        //printf( "%i, %i\n", xi, yi );
        xi = ( xi >= PicSize-1 ) ? PicSize-1 : xi;
        yi = ( yi >= PicSize-1 ) ? PicSize-1 : yi;
        i1 = (int)(( x-h ) / dx);
        i2 = (int)(( x+h ) / dx);
        j1 = (int)(( y-h ) / dy);
        j2 = (int)(( y+h ) / dy);
        if ( i1 != xi || i2 != xi || j1 != yi || j2 != yi ) {
            for ( li=0; li<N; li++ )
                for ( lj=0; lj<N; lj++ ){
                        lx = x + ( li-Nhalf ) * dh;
                        ly = y + ( lj-Nhalf ) * dh;
                        i1 = lx / dx;
                        j1 = ly / dy;
                        if ( i1 < 0 || i1 >= PicSize ||
                                j1 < 0 || j1 >= PicSize ) continue;
                        img[ i1 * PicSize + j1 ] += v * KernelMat2D[0][ li*N + lj ] / ( dx * dy );
                }
        }
        else
            img[ xi * PicSize + yi ] += v / ( dx * dy );
        //printf( "%g\n", img[ xi * PicSize + yi ] );

       image.DataMin = vmin( image.DataMin, v, 1 );
       image.DataMax = vmax( image.DataMax, v );
    }

    for ( i=0; i<SQR(PicSize); i++ ) {
        image.ImgMin = vmin( image.ImgMin, img[i], 1 );
        image.ImgMax = vmax( image.ImgMax, img[i] );
    }

    /*
    MPI_Reduce( &image.DataMin, &image.GlobalDataMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &image.DataMax, &image.GlobalDataMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &image.ImgMin, &image.GlobalImgMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &image.ImgMax, &image.GlobalImgMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );
    */
    find_global_value( image.DataMin, image.GlobalDataMin, MPI_DOUBLE, MPI_MIN );
    find_global_value( image.DataMax, image.GlobalDataMax, MPI_DOUBLE, MPI_MAX );
    find_global_value( image.ImgMin, image.GlobalImgMin, MPI_DOUBLE, MPI_MIN );
    find_global_value( image.ImgMax, image.GlobalImgMax, MPI_DOUBLE, MPI_MAX );

    writelog( "Data:\nMin: %g, Max: %g, GlobalMin: %g, GlobalMax: %g\n"
              "Image: \nMin: %g, Max: %g, GlobalMin: %g, GlobalMax: %g\n",
              image.DataMin, image.DataMax, image.GlobalDataMin, image.GlobalDataMax,
              image.ImgMin,  image.ImgMax,  image.GlobalImgMin,  image.GlobalImgMax );
    image.xmin = All.Start[ proj_i ];
    image.xmax = All.End[ proj_i ];
    image.ymin = All.Start[ proj_j ];
    image.ymax = All.End[ proj_j ];
    image.zmin = All.Start[ proj_k ];
    image.zmax = All.End[ proj_k ];
    writelog( "make slice image ... done\n" );
}
