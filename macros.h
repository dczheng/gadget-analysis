#define writelog( fmt, ... ) { \
    fprintf( LogFileFd, fmt, ##__VA_ARGS__ ); \
    if ( ThisTask == 0 ) { \
        printf( fmt, ##__VA_ARGS__ ); \
    }\
}

#define vmax( a, b ) ( a > b ) ? a : b
#define vmin( a, b, mode) ( mode == 0 ) ? ( ( a > b ) ? b : a ) : ( ( a > b && b > 0 ) ? b : a )
#define check_picture_index( i ) ( ( i<0 || i>=All.PicSize ) ? ( (i<0) ? 0 : All.PicSize ) : i )

#define find_global_value( a, A, type, op ) { \
    MPI_Reduce( &a, &A, 1, type, op, 0, MPI_COMM_WORLD ); \
    MPI_Bcast( &A, 1, type, 0, MPI_COMM_WORLD ); \
    MPI_Barrier( MPI_COMM_WORLD ); \
}

#define check_malloc_var_num() {\
    if ( malloc_n > MALLOC_VAR_NUM ) {\
        writelog( "MALLOC_VAR_NUM IS TOO SMALL ...\n" );\
        endrun( 20180430 ); \
    }\
}

#define check_malloc_var_len( a ) {\
    if ( strlen( #a ) > MALLOC_VAR_LEN ) {\
        writelog( "MALLOC_VAR_LEN IS TOO SMALL ...\n" );\
        endrun( 20180430 ); \
    }\
}

#define malloc_report() { \
    sprintf( malloc_str, "memory info:" );\
    if ( malloc_mem > CUBE( 1024 ) ) {\
        sprintf( malloc_str, "%s Total %g Gb,", malloc_str, malloc_mem / CUBE( 1024. )  );\
    }\
    else if ( malloc_mem > SQR( 1024 ) ) {\
        sprintf( malloc_str, "%s Total %g Mb,", malloc_str, malloc_mem / SQR( 1024. )  );\
    }\
    else if( malloc_mem > 1024 ) {\
        sprintf( malloc_str, "%s Total %g Kb,", malloc_str, malloc_mem / 1024. );\
    }\
    else {\
        sprintf( malloc_str, "%s Total %lli b,", malloc_str, malloc_mem );\
    }\
\
    if ( malloc_max_mem > CUBE( 1024 ) ) {\
        sprintf( malloc_str, "%s Max %g Gb,", malloc_str, malloc_max_mem / CUBE( 1024. ) );\
    }\
    else if ( malloc_max_mem > SQR( 1024 ) ) {\
        sprintf( malloc_str, "%s Max %g Mb,", malloc_str, malloc_max_mem / SQR( 1024. ) );\
    }\
    else if( malloc_max_mem > 1024 ) {\
        sprintf( malloc_str, "%s Max %g Kb,", malloc_str, malloc_max_mem / 1024. );\
    }\
    else {\
        sprintf( malloc_str, "%s Max %lli b,", malloc_str, malloc_max_mem );\
    }\
\
    sprintf( malloc_str, "%s Nvars %lli.\n", malloc_str, malloc_n );\
    writelog( malloc_str ); \
}

#define mymalloc( a, n ) {\
    if ( !(a = malloc( n )) ) { \
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
            writelog( "Failed to allocate memory for `%s` ( %lli b )\n", #a, n ); \
        }\
        endrun( 20180430 ); \
    }\
\
    if ( n > CUBE( 1024 ) ) {\
        writelog( "allocate memory for `%s` ( %g Gb )\n", #a, n / CUBE( 1024. ) ); \
    }\
    else if ( n > SQR( 1024 ) ) {\
        writelog( "allocate memory for `%s` ( %g Mb )\n", #a, n / SQR(  1024. ) ); \
    }\
    else if( n > 1024 ) {\
        writelog( "allocate memory for `%s` ( %g Kb )\n", #a, n /  1024 ); \
    }\
    else {\
        writelog( "allocate memory for `%s` ( %lli b )\n", #a, n ); \
    }\
\
    check_malloc_var_len( a );\
    sprintf( malloc_var[malloc_n], "%s", #a );\
    malloc_var_bytes[malloc_n] = n;\
    malloc_n++; \
    malloc_mem += n; \
    malloc_max_mem = vmax( malloc_max_mem, malloc_mem ); \
    check_malloc_var_num();\
    malloc_report(); \
}

#define myfree( a ) {\
    for ( malloc_i=0; malloc_i<malloc_n; malloc_i++ ) {\
        if ( !( strcmp( malloc_var[malloc_i], #a ) ) ) {\
            malloc_b = malloc_var_bytes[malloc_i];\
            break;\
        }\
    }\
\
    if ( malloc_b > CUBE( 1024 ) ) {\
        writelog( "Free memory for `%s` ( %g Gb )\n", #a, malloc_b / CUBE( 1024. ) ); \
    }\
    else if ( malloc_b > SQR( 1024 ) ) {\
        writelog( "Free memory for `%s` ( %g Mb )\n", #a, malloc_b / SQR(  1024. ) ); \
    }\
    else if( malloc_b > 1024 ) {\
        writelog( "Free memory for `%s` ( %g Kb )\n", #a, malloc_b /  1024 ); \
    }\
    else {\
        writelog( "Free memory for `%s` ( %lli b )\n", #a, malloc_b ); \
    }\
\
    for ( ; malloc_i<malloc_n-1; malloc_i++ ) {\
        sprintf( malloc_var[malloc_i], "%s",  malloc_var[malloc_i+1] ); \
        malloc_var_bytes[malloc_i] = malloc_var_bytes[malloc_i+1]; \
    }\
    malloc_n--; \
    malloc_mem -= malloc_b; \
    malloc_report(); \
}


