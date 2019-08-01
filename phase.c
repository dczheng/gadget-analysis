#include"allvars.h"

void phase() {

    double TempMin, TempMax, DensMin, DensMax,
           LogTempMin, LogTempMax, LogDensMin, LogDensMax,
           DLogTemp, DLogDens, LogDens, LogTemp, Dens, Temp,
           GlobLogTempMin, GlobLogTempMax, GlobLogDensMin, GlobLogDensMax;
    int i, j, k, PicSize;
    char buf[200];

    TempMin = DensMin = DBL_MAX;
    TempMax = DensMax = DBL_MIN;
    writelog( "plot gas phase...\n" );

    //PicSize_tmp = All.PicSize;
    PicSize = All.PicSize;

    reset_img();

    for ( i=0; i<N_Gas; i++ ) {
        vmin20( TempMin, SphP[i].Temp );
        vmax2( TempMax, SphP[i].Temp );
        vmin20( DensMin, SphP[i].Density/All.RhoBaryon );
        vmax2( DensMax, SphP[i].Density/All.RhoBaryon );
    }

    /*
    printf( "%g %g %g %g\n",
        All.PhaseDensMin,
        All.PhaseDensMax,
        All.PhaseTempMin,
        All.PhaseTempMax
        );
    */

    if ( All.PhaseTempMin > 0 )
        TempMin = All.PhaseTempMin;
    if ( All.PhaseTempMax > 0 )
        TempMax = All.PhaseTempMax;

    if ( All.PhaseDensMin > 0 )
        DensMin = All.PhaseDensMin;
    if ( All.PhaseDensMax > 0 )
        DensMax = All.PhaseDensMax;

    LogDensMin = log10( DensMin );
    LogDensMax = log10( DensMax );
    LogTempMin = log10( TempMin );
    LogTempMax = log10( TempMax );
    find_global_value( LogDensMin, GlobLogDensMin, MPI_DOUBLE, MPI_MIN, MpiComm_Master );
    find_global_value( LogDensMax, GlobLogDensMax, MPI_DOUBLE, MPI_MAX, MpiComm_Master );
    find_global_value( LogTempMin, GlobLogTempMin, MPI_DOUBLE, MPI_MIN, MpiComm_Master );
    find_global_value( LogTempMax, GlobLogTempMax, MPI_DOUBLE, MPI_MAX, MpiComm_Master );

    //LogTempMax = 7;

    DLogTemp = ( LogTempMax - LogTempMin ) / (PicSize-1);
    DLogDens = ( LogDensMax - LogDensMin ) / (PicSize-1);

    writelog( "DensMin: %g, DensMax: %g, TempMin: %g, TempMax: %g\n"
            "LogDensMin: %g, LogDensMax: %g, LogTempMin: %g, LogTempMax: %g\n"
            "GlobLogDensMin: %g, GlobLogDensMax: %g,\nGlobLogTempMIn: %g, GlobLogTempMax: %g\n",
            DensMin, DensMax, TempMin, TempMax,
            LogDensMin, LogDensMax, LogTempMin, LogTempMax,
            GlobLogDensMin, GlobLogDensMax, GlobLogTempMin, GlobLogTempMax );

    for ( k=0; k<N_Gas; k++ ) {
        Temp = SphP[k].Temp;
        if ( Temp <= 0 )
            endrun(20181107);
        LogTemp = log10( Temp );

        Dens = SphP[k].Density / All.RhoBaryon;

        if ( Dens <= 0 )
            endrun(20181107);
        LogDens = log10( Dens );

        i = ( LogTemp - LogTempMin ) / DLogTemp;
        check_picture_index( i );

        j = ( LogDens - LogDensMin ) / DLogDens;
        check_picture_index( j );

        //image.img[ i*PicSize + j ] += P[k].Mass / Temp / Dens;
        image.img[ i*PicSize + j ] += P[k].Mass;
  //      image.img[ i*PicSize + j ] ++;
     //  printf( "%g ", image.img[ i*PicSize +j ] ) ;
    }


/*
    for ( i=0, sum=0; i<SQR(PicSize); i++ )
        sum += image.img[i];

    writelog( "sum: %g\n", sum );

    for ( i=0; i<SQR(PicSize); i++ )
        image.img[i] /= sum * DLogTemp * DLogDens;
*/
    for ( i=0; i<SQR(PicSize); i++ )
        image.img[i] /= DLogTemp * DLogDens;

    sprintf( buf, "%sPhase", All.OutputDir );
    create_dir( buf );
    sprintf( buf, "%s/Phase_%03i.dat", buf, All.SnapIndex );

    img_xmin = LogDensMin;
    img_xmax = LogDensMax;
    img_ymin = LogTempMin;
    img_ymax = LogTempMax;
    write_img2( buf, "gas phase" );

    //All.PicSize = PicSize_tmp;
    put_sep;

}
