#include "allvars.h"

void field_to_grid( double *data, double *grid, int pt, int Nmin, int flag ) {

    double fac, v, t, dx[6], *num;
    int i, x[6], NGrid, NGrid2, NGrid3;
    long p, index;

    NGrid = All.NGrid;
    NGrid2 = SQR( NGrid );
    NGrid3 = NGrid2 * NGrid;
    fac = NGrid / All.BoxSize;

    //writelog( "put field to grid ...\n" );
    memset( grid, 0, sizeof(double) * NGrid3 );
    mymalloc2( num, sizeof(double) * NGrid3 );

    for( p=0; p<NumPart6[pt]; p++ ) {
        for( i=0; i<3; i++ ) {

                  t = P[OffsetPart6[pt]+p].Pos[i]* fac;

               x[i] =  t;
             x[i+3] = x[i] + 1;

            dx[i+3] = t - x[i];
            dx[i]   = 1 - dx[i+3];

        }

        for( i=0; i<6; i++ )
            x[i] %= NGrid;

        /*
        index = x[0] * NGrid2 + x[1] * NGrid + x[2];
        grid[index] += data[4*p+3];
        num[index] ++;
        continue;
        */

        //double vv;
        //vv = 0;
        for ( i=0; i<8; i++ ) {

            index = x[(i&1)*3] * NGrid2
                  + x[1+((i>>1)&1)*3] * NGrid 
                  + x[2+((i>>2)&1)*3];

                v = dx[(i&1)*3]
                  * dx[1+((i>>1)&1)*3]
                  * dx[2+((i>>2)&1)*3];
                //vv += v;

            grid[index] += data[p] * v;
            num[index]  += v;

        }
        //printf( "%g\n", vv );

    }

    if ( Nmin > 0 ) 
        for( index=0; index<NGrid3; index++ )
            if ( num[index] < Nmin )
                grid[index] = 0;

    if ( flag )
        for( index=0; index<NGrid3; index++ )
            if ( num[index]>0 )
                grid[index] /= num[index];


    myfree( num );

    //writelog( "put field to grid ... done.\n" );

}
