#define writelog( fmt, ... ) { \
    fprintf( LogFileFd, fmt, ##__VA_ARGS__ ); \
    fflush( LogFileFd ); \
    if ( ThisTask_Master == NTask_Master-1 ) { \
        printf( fmt, ##__VA_ARGS__ ); \
    }\
}

#define vmax( a, b ) ( ( a > b ) ? a : b )
#define vmax2( a, b ) ( a = vmax( a, b ) )
#define vmin( a, b, mode) ( ( mode == 0 ) ? ( ( a > b ) ? b : a ) : ( ( a > b && b > 0 ) ? b : a ) )
#define vmin2( a, b, mode ) ( a = vmin( a, b, mode ) )
#define check_picture_index( i )  i = ( ( i<0 || i>=All.PicSize ) ? ( (i<0) ? 0 : All.PicSize-1 ) : i )
#define PERIODIC( x ) ( ( x > All.HalfBoxSize || x < -All.HalfBoxSize ) ? ( ( x > All.HalfBoxSize ) ? ( x - All.BoxSize ) : ( x + All.BoxSize )  ) : x )
#define NGB_PERIODIC( x ) ( (fabs(x) > All.HalfBoxSize) ? ( All.BoxSize-fabs(x) ) : fabs(x) )

#define DATA_SWAP( a, b, tmp ) {\
   tmp = a; \
   a = b; \
   b = t1; \
}

#define put_sep    writelog( sep_str )
#define put_sep0    writelog( sep_str0 )

#define find_global_value( a, A, type, op ) { \
    MPI_Reduce( &a, &A, 1, type, op, 0, MPI_COMM_WORLD ); \
    MPI_Bcast( &A, 1, type, 0, MPI_COMM_WORLD ); \
    MPI_Barrier( MPI_COMM_WORLD ); \
}

#define check_var_num() {\
    if ( ms.nn > MALLOC_VAR_NUM ) {\
        writelog( "MALLOC_VAR_NUM IS TOO SMALL ...\n" );\
        endrun( 20180430 ); \
    }\
}

#define check_var_len( a ) {\
    if ( strlen( #a ) > MALLOC_VAR_LEN ) {\
        writelog( "MALLOC_VAR_LEN IS TOO SMALL ...\n" );\
        endrun( 20180430 ); \
    }\
}

#define malloc_report() { \
    sprintf( ms.str, "memory info:" );\
    if ( ms.mem > CUBE( 1024 ) ) {\
        sprintf( ms.str, "%s Total %g Gb,", ms.str, ms.mem / CUBE( 1024. )  );\
    }\
    else if ( ms.mem > SQR( 1024 ) ) {\
        sprintf( ms.str, "%s Total %g Mb,", ms.str, ms.mem / SQR( 1024. )  );\
    }\
    else if( ms.mem > 1024 ) {\
        sprintf( ms.str, "%s Total %g Kb,", ms.str, ms.mem / 1024. );\
    }\
    else {\
        sprintf( ms.str, "%s Total %li b,", ms.str, ms.mem );\
    }\
\
    if ( ms.max_mem > CUBE( 1024 ) ) {\
        sprintf( ms.str, "%s Max %g Gb,", ms.str, ms.max_mem / CUBE( 1024. ) );\
    }\
    else if ( ms.max_mem > SQR( 1024 ) ) {\
        sprintf( ms.str, "%s Max %g Mb,", ms.str, ms.max_mem / SQR( 1024. ) );\
    }\
    else if( ms.max_mem > 1024 ) {\
        sprintf( ms.str, "%s Max %g Kb,", ms.str, ms.max_mem / 1024. );\
    }\
    else {\
        sprintf( ms.str, "%s Max %li b,", ms.str, ms.max_mem );\
    }\
\
    sprintf( ms.str, "%s Nvars %li.\n", ms.str, ms.nn );\
    fprintf( MemUseFileFd, ms.str ); \
    fflush( MemUseFileFd ); \
}


#define mymalloc0( a, n, flag, b, shared, disp_unit, mpi_win ) {\
        if ( shared ) \
            MPI_Win_allocate_shared( n, disp_unit, MPI_INFO_NULL,\
                    MpiComm_Local, &a, &mpi_win ); \
        else \
            a = malloc( n ); \
    \
    if ( (!a) && (n > 0) ) { \
        if ( n > CUBE( 1024 ) ) {\
            writelog( "Failed to allocate memory for `%s` ( %g Gb )\n", #a, n / CUBE( 1024. ) ); \
        }\
        else if ( n > SQR( 1024 ) ) {\
            writelog( "Failed to allocate memory for `%s` ( %g Mb )\n", #a, n / SQR( 1024. ) ); \
        }\
        else if( n > 1024 ) {\
            writelog( "Failed to allocate memory for `%s` ( %g Kb )\n", #a, n /  1024. ); \
        }\
        else {\
            writelog( "Failed to allocate memory for `%s` ( %li b )\n", #a, (long)n ); \
        }\
        endrun( 20180430 ); \
    }\
\
    if ( n > CUBE( 1024 ) ) {\
        fprintf( MemUseFileFd, "allocate memory for `%s` ( %g Gb )\n", #a, n / CUBE( 1024. ) ); \
    }\
    else if ( n > SQR( 1024 ) ) {\
        fprintf( MemUseFileFd, "allocate memory for `%s` ( %g Mb )\n", #a, n / SQR(  1024. ) ); \
    }\
    else if( n > 1024 ) {\
        fprintf( MemUseFileFd, "allocate memory for `%s` ( %g Kb )\n", #a, n /  1024. ); \
    }\
    else {\
        fprintf( MemUseFileFd, "allocate memory for `%s` ( %li b )\n", #a, (long)n ); \
    }\
\
    if ( flag == 1 ){\
        fprintf( MemUseFileFd, "initialize `%s` ...\n", #a ); \
        memset( a, b, n ); \
    }\
    check_var_len( a );\
    if ( shared )\
        sprintf( ms.var[ms.nn], "%s", #mpi_win );\
    else \
        sprintf( ms.var[ms.nn], "%s", #a );\
    ms.var_bytes[ms.nn] = n;\
    ms.nn++; \
    ms.mem += n; \
    ms.max_mem = vmax( ms.max_mem, ms.mem ); \
    check_var_num();\
    malloc_report(); \
    fprintf( MemUseFileFd, sep_str ); \
    fflush( MemUseFileFd ); \
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
#define mymalloc2( a, n )     mymalloc( a, n, 1, 0 )
#define mymalloc3( a, n, b )  mymalloc( a, n, 1, b )

#define mymalloc1_shared( a, n, disp_unit, mpi_win )  mymalloc_shared( a, n, 0, 0, disp_unit, mpi_win )
#define mymalloc2_shared( a, n, disp_unit, mpi_win )  mymalloc_shared( a, n, 1, 0, disp_unit, mpi_win )
#define mymalloc3_shared( a, n, b, disp_unit, mpi_win )  mymalloc_shared( a, n, 1, b, disp_unit, mpi_win )

#define myfree0( a, shared, mpi_win ) {\
    for ( ms.i=0; ms.i<ms.nn; ms.i++ ) {\
        if ( shared ) {\
            if ( !( strcmp( ms.var[ms.i], #mpi_win ) ) ) {\
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
    if ( ms.b > CUBE( 1024 ) ) {\
        fprintf( MemUseFileFd, "Free memory for `%s` ( %g Gb )\n", #a, ms.b / CUBE( 1024. ) ); \
    }\
    else if ( ms.b > SQR( 1024 ) ) {\
        fprintf( MemUseFileFd, "Free memory for `%s` ( %g Mb )\n", #a, ms.b / SQR(  1024. ) ); \
    }\
    else if( ms.b > 1024 ) {\
        fprintf( MemUseFileFd, "Free memory for `%s` ( %g Kb )\n", #a, ms.b /  1024. ); \
    }\
    else {\
        fprintf( MemUseFileFd, "Free memory for `%s` ( %li b )\n", #a, ms.b ); \
    }\
\
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
    fprintf( MemUseFileFd, sep_str ); \
    fflush( MemUseFileFd ); \
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
    endrun0( "%s\n", s ); \
}

#define endrun( i ) {\
    endrun0( "%i\n", i ); \
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

