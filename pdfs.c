#include"allvars.h"

void B_Pdf() {

    double BMin, BMax, DensMin, DensMax,
           LogBMin, LogBMax, LogDensMin, LogDensMax,
           DLogB, DLogDens, LogDens, LogB, sum,
           GlobLogBMin, GlobLogBMax, GlobLogDensMin, GlobLogDensMax,
           B;
    int i, j, k, PicSize;
    char buf[200];

    BMin = DensMin = DBL_MAX;
    BMax = DensMax = DBL_MIN;
    writelog( "plot B PDF ...\n" );

    //PicSize_tmp = All.PicSize;
    PicSize = All.PicSize;

    reset_img();
    for ( i=0; i<N_Gas; i++ ) {

        B = get_B( i ) * 1e6;

       // printf( "%g\n", B );
       //
       /*
       if ( B <= 1e-10 )
           continue;
           */

        if ( B == 0 )
            continue;

        vmin20( BMin, B );
        vmax2( BMax, B );

        vmin20( DensMin, SphP[i].Density/All.RhoBaryon);
        vmax2( DensMax, SphP[i].Density/All.RhoBaryon );
    }

    LogDensMin = log10( DensMin );
    LogDensMax = log10( DensMax );
    LogBMin = log10( BMin );
    LogBMax = log10( BMax );
    find_global_value( LogDensMin, GlobLogDensMin, MPI_DOUBLE, MPI_MIN, MpiComm_Master );
    find_global_value( LogDensMax, GlobLogDensMax, MPI_DOUBLE, MPI_MAX, MpiComm_Master );
    find_global_value( LogBMin, GlobLogBMin, MPI_DOUBLE, MPI_MIN, MpiComm_Master );
    find_global_value( LogBMax, GlobLogBMax, MPI_DOUBLE, MPI_MAX, MpiComm_Master );

    //LogBMax = 7;

    DLogB = ( LogBMax - LogBMin ) / (PicSize-1);
    DLogDens = ( LogDensMax - LogDensMin ) / (PicSize-1);
    writelog( "DensMin: %g, DensMax: %g, BMin: %g, BMax: %g\n"
            "LogDensMin: %g, LogDensMax: %g, LogBMin: %g, LogBMax: %g\n"
            "GlobLogDensMin: %g, GlobLogDensMax: %g,\nGlobLogBMIn: %g, GlobLogBMax: %g\n",
            DensMin, DensMax, BMin, BMax,
            LogDensMin, LogDensMax, LogBMin, LogBMax,
            GlobLogDensMin, GlobLogDensMax, GlobLogBMin, GlobLogBMax );

    for ( k=0; k<N_Gas; k++ ) {

        LogB = get_B( k ) * 1e6;

        /*
       if ( B <= 1e-10 )
           continue;
           */
        if ( LogB == 0 )
            continue;

        LogB = log10( LogB );

        //if ( LogB > LogBMax )
        //    continue;

        LogDens = SphP[k].Density / All.RhoBaryon;
        if ( LogDens <= 0 )
            endrun( 20181107 );
        LogDens = log10( LogDens );

        i = ( LogB - LogBMin ) / DLogB;
        check_picture_index( i );

        j = ( LogDens - LogDensMin ) / DLogDens;
        check_picture_index( j );

        image.img[ i*PicSize + j ]++;
     //  printf( "%g ", img[ i*PicSize +j ] ) ;
    }

    for ( i=0, sum=0; i<SQR(PicSize); i++ )
        sum += image.img[i];

    writelog( "sum: %g\n", sum );

    for ( i=0; i<SQR(PicSize); i++ )
        image.img[i] /= sum;

    sprintf( buf, "%sBPdf", All.OutputDir );
    create_dir( buf );
    sprintf( buf, "%s/BPdf_%03i.dat", buf, All.SnapIndex );

    img_xmin = LogDensMin;
    img_xmax = LogDensMax;
    img_ymin = LogBMin;
    img_ymax = LogBMax;
    write_img2( buf, "BPdf" );

    //All.PicSize = PicSize_tmp;
    put_sep;

}

void dens_pdf() {

    double *num, dlogDens, Densmin, Densmax;
    int N, i;
    long p;
    FILE *fd;
    char buf[100];
    
    N = All.DensPdfN;

    get_gas_density_min_max( Densmin, Densmax );
    Densmin /= All.RhoBaryon;
    Densmax /= All.RhoBaryon;

    if ( All.DensPdfMin > 0 )
        Densmin = All.DensPdfMin;
    if ( All.DensPdfMax > 0 )
        Densmax = All.DensPdfMax;

    writelog( "Den, min: %g, max: %g\n", Densmin, Densmax );

    dlogDens = log10( Densmax / Densmin ) / N;
    mymalloc2( num, sizeof(double) * N );

    for( p=0; p<N_Gas; p++ ) {
        i = log10( SphP[p].Density / (Densmin * All.RhoBaryon) ) / dlogDens;
        if ( i > N-1 || i < 0 )
            continue;
        num[i] ++;
    }

    for( i=0; i<N; i++ )
        num[i] /= dlogDens;

    sprintf( buf, "%s/DensPdf", All.OutputDir );
    create_dir( buf );
    sprintf( buf, "%s/DensPdf_%03i.dat", buf, All.SnapIndex );
    fd = fopen( buf, "w" );
    fprintf( fd, "%g %g\n", All.RedShift, All.RedShift );
    for( i=0; i<N; i++ )
        fprintf( fd, "%g %g\n", Densmin*pow( 10, i*dlogDens ), num[i] );

    fclose( fd );
    myfree( num );

}

void T_pdf() {

    double *num, dlogT, Tmin, Tmax;
    int N, i;
    long p;
    FILE *fd;
    char buf[100];
    
    N = All.TPdfN;

    get_gas_temp_min_max( Tmin, Tmax );

    //printf( "%g %g\n", All.TPdfMin, All.TPdfMax );

    if ( All.TPdfMin > 0 )
        Tmin = All.TPdfMin;
    if ( All.TPdfMax > 0 )
        Tmax = All.TPdfMax;
    writelog( "T, min: %g, max: %g\n", Tmin, Tmax );

    dlogT = log10( Tmax / Tmin ) / N;
    mymalloc2( num, sizeof(double) * N );

    for( p=0; p<N_Gas; p++ ) {
        i = log10( SphP[p].Temp / Tmin ) / dlogT;
        if ( i > N-1 || i < 0 )
            continue;
        num[i] ++;
    }

    for( i=0; i<N; i++ )
        num[i] /= dlogT;

    sprintf( buf, "%s/TPdf", All.OutputDir );
    create_dir( buf );
    sprintf( buf, "%s/TPdf_%03i.dat", buf, All.SnapIndex );

    fd = fopen( buf, "w" );
    fprintf( fd, "%g %g\n", All.RedShift, All.RedShift );
    for( i=0; i<N; i++ )
        fprintf( fd, "%g %g\n", Tmin*pow( 10, i*dlogT ), num[i] );

    fclose( fd );
    myfree( num );

}

