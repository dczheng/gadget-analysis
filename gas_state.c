#include"allvars.h"

void gas_state() {

    double TempMin, TempMax, DensMin, DensMax,
           LogTempMin, LogTempMax, LogDensMin, LogDensMax,
           *img, DLogTemp, DLogDens, LogDens, LogTemp, sum,
           GlobLogTempMin, GlobLogTempMax, GlobLogDensMin, GlobLogDensMax;
    int i, j, k, PicSize;
    char buf[200];

    TempMin = DensMin = DBL_MAX;
    TempMax = DensMax = DBL_MIN;
    writelog( "plot gas state...\n" );

    //PicSize_tmp = All.PicSize;
    PicSize = All.PicSize;


    mymalloc2( img, sizeof( double ) * SQR( PicSize ) );

    for ( i=0; i<N_Gas; i++ ) {
        TempMin = vmin( SphP[i].Temp, TempMin, 1 );
        TempMax = vmax( SphP[i].Temp, TempMax );
        DensMin = vmin( SphP[i].Density/All.RhoBaryon, DensMin, 1 );
        DensMax = vmax( SphP[i].Density/All.RhoBaryon, DensMax );
    }

    LogDensMin = log10( DensMin );
    LogDensMax = log10( DensMax );
    LogTempMin = log10( TempMin );
    LogTempMax = log10( TempMax );
    find_global_value( LogDensMin, GlobLogDensMin, MPI_DOUBLE, MPI_MIN );
    find_global_value( LogDensMax, GlobLogDensMax, MPI_DOUBLE, MPI_MAX );
    find_global_value( LogTempMin, GlobLogTempMin, MPI_DOUBLE, MPI_MIN );
    find_global_value( LogTempMax, GlobLogTempMax, MPI_DOUBLE, MPI_MAX );

    //LogTempMax = 7;

    DLogTemp = ( LogTempMax - LogTempMin ) / PicSize;
    DLogDens = ( LogDensMax - LogDensMin ) / PicSize;
    writelog( "DensMin: %g, DensMax: %g, TempMin: %g, TempMax: %g\n"
            "LogDensMin: %g, LogDensMax: %g, LogTempMIn: %g, LogTempMax: %g\n"
            "GlobLogDensMin: %g, GlobLogDensMax: %g,\nGlobLogTempMIn: %g, GlobLogTempMax: %g\n",
            DensMin, DensMax, TempMin, TempMax,
            LogDensMin, LogDensMax, LogTempMin, LogTempMax,
            GlobLogDensMin, GlobLogDensMax, GlobLogTempMin, GlobLogTempMax );

    ZSPRINTF( 0, "1" );
    for ( k=0; k<N_Gas; k++ ) {
        LogTemp = SphP[k].Temp;
        if ( LogTemp < 0 )
            endrun();
        LogTemp = ( LogTemp == 0 ) ? LogTempMin : log10( LogTemp );

        //if ( LogTemp > LogTempMax )
        //    continue;

        LogDens = SphP[k].Density / All.RhoBaryon;
        if ( LogDens < 0 )
            endrun();
        LogDens = ( LogDens == 0 ) ? LogDensMin : log10( LogDens );

        i = ( LogTemp - LogTempMin ) / DLogTemp;
        check_picture_index( i );

        j = ( LogDens - LogDensMin ) / DLogDens;
        check_picture_index( j );
        if ( i < 0 || i >= All.PicSize || j < 0 || j >= All.PicSize)
            printf( "%i, %i\n", i, j );
        img[ i*PicSize + j ]++;
     //  printf( "%g ", img[ i*PicSize +j ] ) ;
    }
    ZSPRINTF( 0, "2" );

    for ( i=0, sum=0; i<SQR(PicSize); i++ )
        sum += img[i];

    for ( i=0; i<SQR(PicSize); i++ )
        img[i] /= sum;

    create_dir( "./gas_state" );
    sprintf( buf, "./gas_state/%s_gas_state_%03i_%.2f.dat", All.FilePrefix, ThisTask, All.RedShift );

    ZSPRINTF( 0, "3" );
    image.img = img;
    img_xmin = LogDensMin;
    img_xmax = LogDensMax;
    img_ymin = LogTempMin;
    img_ymax = LogTempMax;
    write_img2( buf, "gas state" );
    myfree( img );

    //All.PicSize = PicSize_tmp;
    put_block_line;

}
