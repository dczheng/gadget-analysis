#include "allvars.h"

void slice() {
    int pt;
    long offset, num, index, i;
    struct particle_data p_tmp;
    struct sph_particle_data sphp_tmp;
    writelog( "determine slice info ...\n" );

    if ( All.End[All.proj_i] - All.Start[All.proj_i] !=
            All.End[All.proj_j] - All.Start[All.proj_j] ) {
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
}

void make_slice_img( int pt ) {
    double *img, dx, dy, x, y, h, dh, lx, ly, v;
    int i, xi, yi, N, Nhalf, i1, i2, j1, j2, li, lj, PicSize;
    PicSize = All.PicSize;
    writelog( "make slice imgage  ...\n" );
    img = image.img;
    memset( img, 0, sizeof( double ) * PicSize * PicSize );
    dx = dy = (All.End[All.proj_i] - All.Start[All.proj_i])/ PicSize;
    N = All.KernelN;
    Nhalf = N / 2;
    h = All.SofteningTable[pt];
    dh = h / Nhalf;
    for ( i=All.SliceStart[pt]; i<All.SliceEnd[pt]; i++ ){
        x = P[i].Pos[All.proj_i];
        y = P[i].Pos[All.proj_j];
        x -= All.Start[All.proj_i];
        y -= All.Start[All.proj_j];
        v = image.data[ i-All.SliceStart[pt] ];

        //printf( "%g\n", v );
        xi = x / dx;
        yi = y / dy;
        //printf( "%i, %i\n", xi, yi );
        check_picture_index( xi );
        check_picture_index( yi );
        if ( All.KernelInterpolation == 0 ){
            img[ xi * PicSize + yi ] += v;
            continue;
        }
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
                        img[ i1 * PicSize + j1 ] += v * All.KernelMat2D[pt][ li*N + lj ];
                }
        }
        else
            img[ xi * PicSize + yi ] += v;
        //printf( "%g\n", img[ xi * PicSize + yi ] );

    }

    img_xmin = All.Start[ All.proj_i ];
    img_xmax = All.End[ All.proj_i ];
    img_ymin = All.Start[ All.proj_j ];
    img_ymax = All.End[ All.proj_j ];
    writelog( "make slice image ... done\n" );
    put_block_line;

}
