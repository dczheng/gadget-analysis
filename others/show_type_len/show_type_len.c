#include "stdio.h"
#include "limits.h"
#include "float.h"
#include "math.h"

#define ISS( s ) printf( "%-20s: %i  unsigned %-11s: %i\n", #s, sizeof(s), #s, sizeof(unsigned s) );
#define ISS2( s1, s2 ) printf( "%-20s: %i  unsigned %-11s: %i\n", #s1" "#s2, sizeof(s1 s2), #s1" "#s2, sizeof(unsigned s1 s2) );
#define FSS( s ) printf( "%-20s: %i\n", #s, sizeof(s) );
#define FSS2( s1, s2 ) printf( "%-20s: %i\n", #s1" "#s2, sizeof(s1 s2) );

void main() {
    printf("\n");
    ISS( char );
    ISS( int );
    ISS( short );
    ISS( long );
    ISS2( long, long );
    FSS( float );
    FSS( double );
    FSS2( long, double );
    printf( "\n" );

    printf("CHAR_BIT   = %d\n", CHAR_BIT);
    printf("MB_LEN_MAX = %d\n", MB_LEN_MAX);
    printf("\n");

    printf("CHAR_MIN   = %+d\n", CHAR_MIN);
    printf("CHAR_MAX   = %+d\n", CHAR_MAX);
    printf("SCHAR_MIN  = %+d\n", SCHAR_MIN);
    printf("SCHAR_MAX  = %+d\n", SCHAR_MAX);
    printf("UCHAR_MAX  = %u\n",  UCHAR_MAX);
    printf("\n");

    printf("SHRT_MIN   = %+d\n", SHRT_MIN);
    printf("SHRT_MAX   = %+d\n", SHRT_MAX);
    printf("USHRT_MAX  = %u\n",  USHRT_MAX);
    printf("\n");

    printf("INT_MIN    = %+d\n", INT_MIN);
    printf("INT_MAX    = %+d\n", INT_MAX);
    printf("UINT_MAX   = %u\n",  UINT_MAX);
    printf("\n");

    printf("LONG_MIN   = %+ld\n", LONG_MIN);
    printf("LONG_MAX   = %+ld\n", LONG_MAX);
    printf("ULONG_MAX  = %lu\n",  ULONG_MAX);
    printf("\n");

    printf("LLONG_MIN  = %+lld\n", LLONG_MIN);
    printf("LLONG_MAX  = %+lld\n", LLONG_MAX);
    printf("ULLONG_MAX = %llu\n",  ULLONG_MAX);
    printf("\n");

    printf("FLT_RADIX    = %d\n", FLT_RADIX);
    printf("FLT_MIN      = %e\n", FLT_MIN);
    printf("FLT_MAX      = %e\n", FLT_MAX);
    printf("FLT_EPSILON  = %e\n", FLT_EPSILON);
    printf("FLT_DIG      = %d\n", FLT_DIG);
    printf("FLT_MIN_EXP  = %d\n",  FLT_MIN_EXP);
    printf("FLT_MIN_10_EXP  = %d\n",  FLT_MIN_10_EXP);
    printf("FLT_MAX_EXP     = %d\n",  FLT_MAX_EXP);
    printf("FLT_MAX_10_EXP  = %d\n",  FLT_MAX_10_EXP);
    //printf("FLT_ROUNDS      = %d\n",  FLT_ROUNDS);
    printf( "\n" );

    printf("DBL_MIN      = %e\n", DBL_MIN);
    printf("DBL_MAX      = %e\n", DBL_MAX);
    printf("DBL_EPSILON  = %e\n", DBL_EPSILON);
    printf("DBL_DIG      = %d\n", DBL_DIG);
    printf("DBL_MIN_EXP  = %d\n",  DBL_MIN_EXP);
    printf("DBL_MIN_10_EXP  = %d\n",  DBL_MIN_10_EXP);
    printf("DBL_MAX_EXP     = %d\n",  DBL_MAX_EXP);
    printf("DBL_MAX_10_EXP  = %d\n",  DBL_MAX_10_EXP);
    printf( "\n" );

    printf("LDBL_MIN      = %Le\n", LDBL_MIN);
    printf("LDBL_MAX      = %Le\n", LDBL_MAX);
    printf("LDBL_EPSILON  = %Le\n", LDBL_EPSILON);
    printf("LDBL_DIG      = %d\n", LDBL_DIG);
    printf("LDBL_MIN_EXP  = %d\n",  LDBL_MIN_EXP);
    printf("LDBL_MIN_10_EXP  = %d\n",  LDBL_MIN_10_EXP);
    printf("LDBL_MAX_EXP     = %d\n",  LDBL_MAX_EXP);
    printf("LDBL_MAX_10_EXP  = %d\n",  LDBL_MAX_10_EXP);

}
