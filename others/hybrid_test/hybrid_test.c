#include "omp.h"
#include "stdio.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "time.h"
#include "gsl/gsl_integration.h"
#include "sys/stat.h"
#include "limits.h"

char *sep_str = "---------------------------------------\n";
char *sep_str2 = "#######\n";
#define put_sep printf( sep_str )
#define put_sep2 printf( sep_str2 )

//#define HYBRID

#ifdef HYBRID
int ThisTask, NTask, provided;
#endif

double trapzd( double (*func)( double, void* ),
        void *params, double a, double b, int n ) {

    int i, num;
    double sum, step, x;
    static double s;

    if ( n == 1 ) {
        s = 0.5 * ( (*func)( a, params ) + (*func)(b, params) ) * (b-a);
        return s;
    }

    for( num=1, i=1; i<n-1; i++ ) num <<= 1;

    step = (b-a) / num;

    x = a+0.5*step;

    for( sum=0, i=0; i<num; i++, x+=step )
        sum += (*func)(x, params);

    s = 0.5 * ( s + (b-a)*sum / num );

    return s;

}

#define REFINE_MAX 20
double qtrap( double (*func)( double, void* ),
        void *params, double a, double b, double err ) {

    double s, olds;
    int i;

    olds = -DBL_MAX;

    for( i=1; i<=REFINE_MAX; i++ ) {

        s = trapzd( func, params, a, b, i );

        if ( i>5 )
            if ( fabs( s-olds ) < err*fabs(olds) ||
                    ( s == 0 && olds == 0 ) )
                return s;

        olds = s;
    }

    printf( "a: %g, b: %g\n", a, b );
    printf( "Too many steps in qtrap!\n" );
    exit( 0 );

}

double second() {
    return ( (double)clock() / CLOCKS_PER_SEC );
}

double myF( double x, void *params ) {
    return x*x;
}

void my_compute() {

    int i;

    for( i=1; i<1e3; i++ )
        qtrap( &myF, NULL, 1, 10, 1e-1 );

    /*
    double r;
    (void)r;
    for( i=1; i<1e7; i++ )
        //myF( 1.1, NULL );
        r = 1.1 * 1.1;
        */


    /*
    double r, err;
    gsl_integration_workspace *inte_ws =
        gsl_integration_workspace_alloc( 1000 );
    gsl_function F;
    F.function = &_F;
    F.params = NULL;

    gsl_integration_qag( &F, 1, 100, 0, 1e-1, 1000, GSL_INTEG_GAUSS61,
            inte_ws, &r, &err );

    gsl_integration_workspace_free( inte_ws );
    */

}


void omp_test0() {

    int i, num;

#pragma omp parallel private(i, num)
{

    num = 10000;
    int thread_id, thread_num;
    thread_id = omp_get_thread_num();
    thread_num = omp_get_num_threads();

    for ( i=0; i<num; i++ ) {

       /*
        double t;
        t = omp_get_wtime();
        */

        if ( i % thread_num ==  thread_id )
            my_compute();

       /*
       printf( "[%2i], %i, %g\n", omp_get_thread_num(), i,
               omp_get_wtime() - t );
               */

    }
}

}

#ifdef HYBRID
void hybrid_test( int argc, char **argv ) {

    printf( "hybrid test ...\n" );
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &provided );
    //MPI_Init_thread( &argc, &argv, MPI_THREAD_SERIALIZED, &provided );
    MPI_Comm_rank( MPI_COMM_WORLD, &ThisTask );
    MPI_Comm_size( MPI_COMM_WORLD, &NTask );

    MPI_Barrier( MPI_COMM_WORLD );

    MPI_Finalize();

}
#endif


void omp_test() {

    int num;
    double timer1;

    put_sep;
    printf( "omp test ...\n" );
    put_sep;
    timer1 = omp_get_wtime();

    num = 1;
    omp_set_num_threads( num );

    omp_test0();
    printf( "mean time[%i]: %g\n", num, (omp_get_wtime() - timer1) );
    put_sep;

    num = 5;
    omp_set_num_threads( num );

    timer1 = omp_get_wtime();

    omp_test0();
    printf( "mean time[%i]: %g\n", num, (omp_get_wtime() - timer1) );
    put_sep;

}

int main( int argc, char **argv ) {

#ifdef HYBRID
    (void) ThisTask;
    (void) NTask;
    (void) provided;
    hybrid_test( argc, argv );
#else
    omp_test();
#endif


    return 0;

}
