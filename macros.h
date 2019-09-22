#define writelog( fmt, ... ) { \
    fprintf( LogFileFd, fmt, ##__VA_ARGS__ ); \
    fflush( LogFileFd ); \
    if ( ThisTask_Master == NTask_Master-1 ) { \
        printf( fmt, ##__VA_ARGS__ ); \
    }\
}

#define myfopen( opt, fmt, ... ) ({\
    char myfopenbuf[100];\
    FILE *fd;\
    do {\
        sprintf( myfopenbuf, fmt, ##__VA_ARGS__ );\
        fd = fopen( myfopenbuf, opt );\
        if ( NULL == fd ) {\
            printf( "can not open `%s`.\n", myfopenbuf );\
            endrun(20190826);\
        }\
    }while(0);\
    fd;\
})

#define create_dir( fmt, ... ) {\
    char cbuf[100];\
    sprintf( cbuf, fmt, ##__VA_ARGS__ );\
    create_dir0( cbuf );\
}

#define vmax( a, b ) ( ( (a) > (b) ) ? (a) : (b) )
#define vmax2( a, b ) ( a = vmax( a, b ) )
#define vmin( a, b, mode) ( ( mode == 0 ) ? ( ( (a) > (b) ) ? (b) : (a) ) : ( ( (a) > (b) && (b) > 0 ) ? (b) : (a) ) )
#define vmin2( a, b ) ( a = vmin( a, b, 0 ) )
#define vmin20( a, b ) ( a = vmin( a, b, 1 ) )
#define check_picture_index( i )  i = ( ( (i)<0 || (i)>=PicSize ) ? ( ((i)<0) ? 0 : PicSize-1 ) : (i) )

#define PERIODIC_HALF( x ) ( ( (x) > HalfBoxSize || (x) < -HalfBoxSize ) ? ( ( (x) > HalfBoxSize ) ? ( (x) - BoxSize ) : ( x + BoxSize )  ) : (x) )
#define PERIODIC( x ) ( ( (x) > BoxSize || (x) <  0) ? ( ( (x) > BoxSize ) ? ( (x) - BoxSize ) : ( x + BoxSize )  ) : (x) )
#define NGB_PERIODIC( x ) ( (fabs(x) > HalfBoxSize) ? ( BoxSize-fabs(x) ) : fabs(x) )

#define put_sep    writelog( sep_str )
#define put_sep0    writelog( sep_str0 )

#define check_var_num() {\
    if ( ms.nn > MALLOC_VAR_NUM ) {\
        endruns( "MALLOC_VAR_NUM IS TOO SMALL ..." );\
    }\
}

#define check_var_len( a ) {\
    if ( strlen( #a ) > MALLOC_VAR_LEN ) {\
        endruns( "MALLOC_VAR_LEN IS TOO SMALL ...\n" );\
    }\
}

#define fmt_mem_used( mm ) {\
    if ( mm > CUBE( 1024 ) ) {\
        sprintf( ms.mem_str, "%g Gb", mm / CUBE( 1024. )  );\
    }\
    else if ( mm > SQR( 1024 ) ) {\
        sprintf( ms.mem_str, "%g Mb", mm / SQR( 1024. )  );\
    }\
    else if( mm > 1024 ) {\
        sprintf( ms.mem_str, "%g Kb", mm / 1024. );\
    }\
    else {\
        sprintf( ms.mem_str, "%li b", (long)mm );\
    }\
\
}

#define malloc_report() { \
    fprintf( UsedMemFileFd, "Memory info:\n" ); \
    fprintf( UsedMemFileFd, sep_str0 ); \
    for ( ms.i=0; ms.i<ms.nn; ms.i++ ) {\
        fmt_mem_used( ms.var_bytes[ms.i] );\
        fprintf( UsedMemFileFd, "%-20s: %s\n", ms.var[ms.i], ms.mem_str ); \
    }\
    fprintf( UsedMemFileFd, sep_str0 ); \
    fmt_mem_used( ms.mem );\
    sprintf( ms.str, "Total %s,", ms.mem_str  );\
    fmt_mem_used( ms.max_mem );\
    sprintf( ms.str, "%s Max %s,", ms.str, ms.mem_str  );\
    sprintf( ms.str, "%s Nvars %li.\n", ms.str, ms.nn );\
    fprintf( UsedMemFileFd, ms.str ); \
    fprintf( UsedMemFileFd, sep_str ); \
    fflush( UsedMemFileFd ); \
}


#define mymalloc0( a, n, flag, b, shared, disp_unit, mpi_win ) {\
    fprintf( UsedMemFileFd, "%s %s %i\n", __FILE__, __FUNCTION__, __LINE__ ); \
        if ( shared ) \
            MPI_Win_allocate_shared( n, disp_unit, MPI_INFO_NULL,\
                    MpiComm_Local, &a, &mpi_win ); \
        else \
            a = malloc( n ); \
    \
    fmt_mem_used( n );\
    if ( (!a) && (n > 0) ) { \
        endrun0( "Failed to allocate memory for `%s` ( %s )\n", #a, ms.mem_str ); \
    }\
    fprintf( UsedMemFileFd, "allocate memory for `%s` ( %s )\n", #a, ms.mem_str ); \
    if ( flag == 1 ){\
        fprintf( UsedMemFileFd, "initialize `%s` ...\n", #a ); \
        memset( a, b, n ); \
    }\
    check_var_len( a );\
    if ( shared )\
        sprintf( ms.var[ms.nn], "%s", &((#mpi_win)[7]) );\
    else \
        sprintf( ms.var[ms.nn], "%s", #a );\
    ms.var_bytes[ms.nn] = n;\
    ms.nn++; \
    ms.mem += n; \
    ms.max_mem = vmax( ms.max_mem, ms.mem ); \
    check_var_num();\
    malloc_report(); \
}

#define mymalloc( a, n, flag, b ) mymalloc0( a, n, flag, b, 0, ThisTask, MpiWin_P ) // here ThisTask and MpiWin_P are useless.
#define mymalloc_shared( a, n, flag, b, disp_unit, mpi_win ) { \
    if ( ThisTask_Local == 0 ) { \
        mymalloc0( a, n, flag, b, 1, disp_unit, mpi_win ) \
    } \
    else { \
        mymalloc0( a, 0, flag, b, 1, disp_unit, mpi_win ) \
        MPI_Win_shared_query( mpi_win, 0, &WinSize, &WinDisp, &a ); \
    }\
}

#define mymalloc1( a, n )     mymalloc( a, n, 0, 0 )
#define mymalloc2( a, n )     mymalloc( a, n, 1, 0 )  // initial to 0
#define mymalloc3( a, n, b )  mymalloc( a, n, 1, b )

#define mymalloc1_shared( a, n, disp_unit, mpi_win )  mymalloc_shared( a, n, 0, 0, disp_unit, mpi_win )
#define mymalloc2_shared( a, n, disp_unit, mpi_win )  mymalloc_shared( a, n, 1, 0, disp_unit, mpi_win )
#define mymalloc3_shared( a, n, b, disp_unit, mpi_win )  mymalloc_shared( a, n, 1, b, disp_unit, mpi_win )

#define myfree0( a, shared, mpi_win ) {\
    fprintf( UsedMemFileFd, "%s %s %i\n", __FILE__, __FUNCTION__, __LINE__ ); \
    for ( ms.i=ms.nn-1; ms.i>=0; ms.i-- ) {\
        if ( shared ) {\
            if ( !( strcmp( ms.var[ms.i], &((#mpi_win)[7]) ) ) ) {\
                ms.b = ms.var_bytes[ms.i];\
                break;\
            }\
        }\
        else {\
            if ( !( strcmp( ms.var[ms.i], #a ) ) ) {\
                ms.b = ms.var_bytes[ms.i];\
                break;\
            }\
        }\
    }\
\
    fmt_mem_used( ms.b );\
    fprintf( UsedMemFileFd, "Free memory for `%s` ( %s )\n", ms.var[ms.i], ms.mem_str ); \
    for ( ; ms.i<ms.nn-1; ms.i++ ) {\
        sprintf( ms.var[ms.i], "%s",  ms.var[ms.i+1] ); \
        ms.var_bytes[ms.i] = ms.var_bytes[ms.i+1]; \
    }\
    ms.nn--; \
    ms.mem -= ms.b; \
\
    if ( shared ) \
        MPI_Win_free( &mpi_win );\
    else \
        free( a ); \
\
    malloc_report(); \
}

#define myfree( a )        myfree0( a, 0, MpiWin_P )
#define myfree_shared( mpi_win ) myfree0( P, 1, mpi_win )

#define endrun0( fmt, ... ) {\
    fprintf( stderr, fmt, ##__VA_ARGS__ ); \
    fprintf( stderr, "END IN: ( %s %s %i )\n" , __FILE__, __FUNCTION__, __LINE__ ); \
    MPI_Abort( MPI_COMM_WORLD, 0 ); \
    exit( 0 ); \
}

#define endruns( s ) {\
    endrun0( "error info: %s\n", s ); \
}

#define endrun( i ) {\
    endrun0( "error level: %i\n", i ); \
}


#define mytimer_start() \
    double timer1, timer2;\
    (void) timer1;\
    (void) timer2;\
writelog( "[Timer Start in `%s`]\n", __FUNCTION__ ); \
    timer1 = second(); \
    timer2 = timer1;

#define mytimer() \
    writelog( "[Time in `%s`]: %g sec\n", __FUNCTION__, second() - timer2 ); \
    timer2 = second();

#define mytimer_end() \
    writelog( "[Total Time in `%s`]: [%g sec]\n", __FUNCTION__, second() - timer1 ); \

#define PDF2D_BIT_MODE          1 
#define PDF2D_BIT_XLOG          2 
#define PDF2D_BIT_YLOG          4 
#define PDF2D_BIT_FIXEDX        8 
#define PDF2D_BIT_FIXEDY        16 
#define PDF2D_BIT_NORM          32 
#define PDF2D_BIT_UNIT_AREA     64 

#define pdf2d( x, y, w, num, dn, flag, mm ) {\
    if ( flag & PDF2D_BIT_MODE ) \
        flag -= PDF2D_BIT_MODE;\
    pdf2d_or_field2d( x, y, w, num, dn, flag, mm, 0 );\
}

#define field2d( x, y, z, num, dn, flag, mm, Nmin ) {\
    if ( !(flag & PDF2D_BIT_MODE) ) \
        flag += PDF2D_BIT_MODE;\
    pdf2d_or_field2d( x, y, z, num, dn, flag, mm, Nmin );\
}

#define get_index( i, j, k, N ) ( ((i)*(N)+(j))*(N)+(k) )
#define get_B( i ) ( sqrt(SQR(SphP[i].B[0]) + SQR(SphP[i].B[1]) + SQR(SphP[i].B[2]) ) )

#define get_pressure( i ) (SphP[i].u*( GAMMA_MINUS1*SphP[i].Density/Time3 ))

#define get_gas_density_min_max( dmin, dmax ) {\
    long p;\
    dmin = DBL_MAX;\
    dmax = -dmin;\
    for( p=0; p<N_Gas; p++ ) {\
        vmax2( dmax, SphP[p].Density );\
        vmin2( dmin, SphP[p].Density );\
        }\
}

#define get_B_min_max( bmin, bmax ) {\
    long p;\
    double b;\
    bmin = DBL_MAX;\
    bmax = -bmin;\
    for( p=0; p<N_Gas; p++ ) {\
        b = get_B( p );\
        vmax2( bmax, b );\
        vmin20( bmin, b );\
        }\
}

#define get_gas_temp_min_max( Tmin, Tmax ) {\
    long p;\
    Tmin = DBL_MAX;\
    Tmax = -Tmin;\
    for( p=0; p<N_Gas; p++ ) {\
        vmax2( Tmax, SphP[p].Temp );\
        vmin2( Tmin, SphP[p].Temp );\
        }\
}

#define make_group_output_filename( buf, nstr, group_index ) \
    sprintf( buf, "%s%s/%s_%03i_%04i_%c.dat",\
            GroupDir, nstr, nstr, SnapIndex, group_index, Sproj );

#define put_header( s ) {\
    writelog( ">>> %s\n", s );\
    WATCH_POINT( "debug point" );\
}
#define put_end() {\
    WATCH_POINT( "debug point" );\
}

#define CRE_F( c, a, q1, q2, q ) ( ((q) < (q1) || (q) > (q2) ) ? 0 : (c) * pow((q), -(a)) )
