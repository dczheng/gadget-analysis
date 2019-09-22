
#ifdef ZDEBUG

void signal_hander( int sig );
void init_zsig();
void free_zsig();
void empty_buf_zsig();
void stack_push_zsig( double a );
void stack_pop_zsig();
void stack_reset_zsig();

MyIDType get_id( int i );

#include "signal.h"
#include "string.h"
#define ZDEBUG_NUM 10
#define ZDEBUG_BUF_LEN 30720

struct zsig_struct{
    char stop[ZDEBUG_BUF_LEN];
    char buf[ZDEBUG_NUM][ZDEBUG_BUF_LEN];
    char str_fmt[ZDEBUG_BUF_LEN];
    long i, j, stack_max, stack_pos, no;
    double *stack;
} zsig;

#define RAISE_SIGSTOP() { \
    sprintf( zsig.stop, "%s, %s, %i" , __FILE__, __FUNCTION__, __LINE__ ); raise( SIGSEGV );\
}

#define ZSPRINTF( k, fmt, ... ) { \
    if ( k >= ZDEBUG_NUM ) { \
        sprintf( zsig.stop, "buffer number is too few." ); \
        RAISE_SIGSTOP(); \
    } \
    sprintf( zsig.str_fmt, "%s #( %s, %s, %i )", fmt, __FILE__, __FUNCTION__, __LINE__ ); \
    sprintf( zsig.buf[k], zsig.str_fmt, ##__VA_ARGS__ ); \
}

#define ZSPRINTFA( k, fmt, ... ) { \
    if ( k >= ZDEBUG_NUM ) { \
        sprintf( zsig.stop, "buffer number is too few." ); \
        RAISE_SIGSTOP(); \
    } \
    if ( strlen( zsig.buf[k] ) > ZDEBUG_BUF_LEN ) { \
        sprintf( zsig.stop, "buffer is too small." ); \
        RAISE_SIGSTOP(); \
    }\
    sprintf( zsig.str_fmt, "%s\n%s #( %s, %s, %i )", zsig.buf[k], fmt, __FILE__, __FUNCTION__, __LINE__ ); \
    sprintf( zsig.buf[k], zsig.str_fmt, ##__VA_ARGS__ ); \
}

#define WATCH_POINT( ... ) { \
    sprintf( zsig.stop, "%s, %s, %i: %s", __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__ ); \
}

#define CHECK_NAN3( A ) { \
    for ( zsig.i=0; zsig.i<3; zsig.i++ ) { \
        if ( isnan( A[ zsig.i ] ) ) \
            RAISE_SIGSTOP(); \
    } \
}

#define CHECK_NAN( A ) { \
    if ( isnan( A ) ) \
            RAISE_SIGSTOP(); \
}

#define CHECK_INF( A ) { \
    if ( isinf( A ) ) \
            RAISE_SIGSTOP(); \
}

#define CHECK_NAN_INF( A ) {\
    CHECK_NAN( A ); \
    CHECK_INF( A ); \
}
#define CHECK_NAN_INF_NEG( A ) {\
    CHECK_NAN( A ); \
    CHECK_INF( A ); \
    CHECK_NEG( A ); \
}

#define CHECK_NEG( A ) { \
    if ( A < 0 ) \
            RAISE_SIGSTOP(); \
}

#define CHECK_LARGE_NUMBER( A, B ) { \
    if ( A >= B ) \
            RAISE_SIGSTOP(); \
}

#define DEBUG_CHECK_TEMP() {\
    double fac_entr_to_temp, T;\
    for( zsig.i=0; zsig.i<N_gas; zsig.i++ ) { \
    fac_entr_to_temp = pow( SphP[zsig.i].d.Density  / (All.Time * All.Time * All.Time), GAMMA_MINUS1 ) / GAMMA_MINUS1 \
    *  (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN     * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g; \
        T = SphP[zsig.i].Entropy * fac_entr_to_temp;\
        if ( T > 1e7 ) {\
            ZSPRINTF( 0, "T: %g, Entropy: %g, EntropyPred: %g, DtEntropy: %g\n",\
            T, SphP[zsig.i].Entropy, SphP[zsig.i].EntropyPred, SphP[zsig.i].e.DtEntropy);\
            RAISE_SIGSTOP();\
        } \
    }\
}

#define DEBUG_CHECK_TEMP_SINGLE( ii ) {\
    double fac_entr_to_temp, T, ha, RhoBar;\
    ha = hubble_function(All.Time);\
    RhoBar = All.OmegaBaryon * 3 * ha * ha / ( 8*M_PI*All.G );\
    fac_entr_to_temp = pow( SphP[ii].d.Density  / (All.Time * All.Time * All.Time), GAMMA_MINUS1 ) / GAMMA_MINUS1 \
    *  (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN     * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g; \
        T = SphP[ii].Entropy * fac_entr_to_temp;\
        if ( T > 1e7 ) {\
            ZSPRINTF( 0, "T: %g, Entropy: %g, EntropyPred: %g, DtEntropy: %g, Dens: %g, ngb: %f\n", T, SphP[ii].Entropy, SphP[ii].EntropyPred, SphP[ii].e.DtEntropy, SphP[ii].d.Density/RhoBar, SphP[ii].n.NumNgb);\
            printf( "T: %g, Entropy: %g, EntropyPred: %g, DtEntropy: %g, Dens: %g, ngb: %f\n", T, SphP[ii].Entropy, SphP[ii].EntropyPred, SphP[ii].e.DtEntropy, SphP[ii].d.Density/RhoBar, SphP[ii].n.NumNgb);\
        } \
}

#define ZCLEAR()   { empty_buf_zsig();  }
#define ZPUSH( a ) { stack_push_zsig( (double)a ); }
#define ZPOP()     { stack_pop_zsig(); }
#define ZRESET()   { stack_reset_zsig(); }
#define CHECK_TREE( a, lim ) { \
    ZRESET(); \
    ZPUSH( All.MaxPart ); \
    zsig.no = a;\
    do {\
        ZPUSH( zsig.no ); \
        zsig.no = Nodes[zsig.no].u.d.father; \
        if ( zsig.no > lim ) {\
            ZPUSH( zsig.no ); \
            RAISE_SIGSTOP(); \
        }\
    }while( zsig.no>=0 ); \
}

#define CHECK_FATHER( no, a, mode ) { \
    if ( mode > 0 ) \
        if ( Nodes[no].u.d.father > a )\
            RAISE_SIGSTOP(); \
    if ( mode == 0 ) \
        if ( Nodes[no].u.d.father == a )\
            RAISE_SIGSTOP(); \
    if ( mode < 0 ) \
        if ( Nodes[no].u.d.father < a )\
            RAISE_SIGSTOP(); \
}

#else

#define RAISE_SIGSTOP()
#define ZSPRINTF( k, fmt, ... )
#define ZSPRINTFA( k, fmt, ... )
#define WATCH_POINT( ... )
#define CHECK_NAN3( A )
#define CHECK_NAN( A )
#define CHECK_NEG( A )
#define CHECK_LARGE_NUMBER( A, B )
#define ZCLEAR()
#define ZPUSH( a )
#define ZPOP()
#define ZRESET()
#define CHECK_TREE( no, lim )
#define CHECK_FATHER( no, a )
#define CHECK_INF( A )
#define DEBUG_CHECK_TEMP()
#define DEBUG_CHECK_TEMP_SINGLE( ii )
#define CHECK_NAN_INF( A )

#endif
