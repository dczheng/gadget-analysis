#include "allvars.h"

double kernel( double q ) {
    if ( 0 <= q && q <= 0.5 )
        return 8.0 / PI * ( 1 - 6*q*q + 6*q*q*q );
    if ( q <= 1 )
        return 8.0 / PI * ( 2.0 * pow( 1-q, 3.0 ) );
    return 0;
}

void init_kernel_matrix() {
    int pt, i, j, k, N, Nhalf;
    long num;
    double dh, q, h;
    writelog( "initialize kernel matrix ...\n" );
    N = All.KernelN;
    Nhalf = N / 2;
    for( pt=0; pt<6; pt++ ) {
        num = header.npartTotal[pt] + ( ( ( long )header.npartTotalHighWord[pt] ) << 32 );
        if ( num != 0 ) {
            KernelMat2D[pt] = malloc( sizeof( double ) * N * N );
            KernelMat3D[pt] = malloc( sizeof( double ) * N * N * N );
            memset( KernelMat2D[pt], 0, sizeof( double ) * N * N );
            memset( KernelMat3D[pt], 0, sizeof( double ) * N * N * N );

            if ( SofteningTable[pt] == 0 )
                endrun0( "SofteningTable[%d] is zeros !!!\n", pt );

            h = SofteningTable[pt];
            dh = h / Nhalf;
            for ( i=0; i<N; i++ )
                for ( j=0; j<N; j++ )
                    for ( k=0; k<N; k++ ) {
                        //writelog( "%i\n", Nhalf );
                        q = sqrt( SQR( i-Nhalf ) +
                                  SQR( j-Nhalf ) +
                                  SQR( k-Nhalf ) );
                        q = q / Nhalf;
                        KernelMat3D[pt][ i*N*N + j*N + k ] = kernel( q ) / pow( h, 3 )
                            * pow( dh, 3 );
                        KernelMat2D[pt][ i*N + j ] += KernelMat3D[pt][ i*N*N + j*N + k ];
                        //writelog( "%f, %f\n", q, KernelMat2D[pt][i*N+j] );
                    }
            /*
            double t;
            t = 0;
            for ( i=0; i<N; i++ )
                for ( j=0; j<N; j++ ) {
                    t += KernelMat2D[pt][ i*N + j ];
                }
            writelog( "%f\n", t );
            endrun( 0 );
            */

            for ( i=0, q=0; i<N*N*N; i++ ) {
                q += KernelMat3D[pt][i];
            }
            writelog( "sum of 3d kernel[%d] = %g\n", pt, q );
            for ( i=0, q=0; i<N*N; i++ )
                q += KernelMat2D[pt][i];
            writelog( "sum of 2d kernel[%d] = %g\n", pt, q );
        }
    }
    writelog( "initialize kernel matrix ... done.\n" );
}

void free_kernel_matrix() {
    long num;
    int pt;
    writelog( "free kernel matrix ...\n" );
    for ( pt=0; pt<6; pt++ ) {
        num = header.npartTotal[pt] + ( ( ( long )header.npartTotalHighWord[pt] ) << 32 );
        if ( num != 0 ) {
            free( KernelMat2D[pt] );
            free( KernelMat3D[pt] );
        }
    }
}

#ifndef KERNEL_H
#define KERNEL_H

/* fall back to standard kernel */
#if !defined(QUINTIC_KERNEL) && !defined(WENDLAND_C4_KERNEL) && !defined(WENDLAND_C6_KERNEL)
#define STANDARD_KERNEL   
#endif

/* fall back to three dimensions */
#if !defined(TWODIMS) && !defined(ONEDIM)
#define THREEDIMS
#endif

/* Norms */
#ifdef STANDARD_KERNEL
#ifdef THREEDIMS
#define  NORM 8.0/M_PI                   /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define  NORM 40.0/(7.0*M_PI)	         /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIM
#define  NORM  4.0/3.0        	         /*!< For 1D-normalized kernel */
#endif
#endif /* STANDARD_KERNEL */

#ifdef QUINTIC_KERNEL
#ifdef THREEDIMS
#define  NORM 2187.0/(40.0*M_PI)	    /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define  NORM 15309.0/(478.0*M_PI)	    /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIM
#define  NORM 243.0/40.0        	    /*!< For 1D-normalized kernel */
#endif
#endif  /* QUINTIC_KERNEL */

#ifdef WENDLAND_C4_KERNEL
#ifdef THREEDIMS
#define  NORM 495.0/(32*M_PI)	        /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define  NORM 9.0/M_PI	                /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIM
#define  NORM 3.0/2.0      	            /*!< For 1D-normalized kernel */
#endif
#endif  /* WENDLAND_C4_KERNEL */

#ifdef WENDLAND_C6_KERNEL
#ifdef THREEDIMS
#define  NORM 1365.0/(64*M_PI)	        /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define  NORM 78.0/(7*M_PI)	                /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIM
#define  NORM 55.0/32.0      	            /*!< For 1D-normalized kernel */
#endif
#endif  /* WENDLAND_C6_KERNEL */

inline void kernel_hinv(double h, double *hinv, double *hinv3, double *hinv4)
{
  *hinv = 1.0 / h;

#ifdef THREEDIMS
  *hinv3 = *hinv * *hinv * *hinv;
#endif

#ifdef TWODIMS
  *hinv3 = *hinv * *hinv / boxSize_Z; /* typo ?? */
#endif

#ifdef ONEDIM
  *hinv3 = *hinv;
#endif

  *hinv4 = *hinv3 * *hinv;

  return;
} 

/* Attention: Here we assume that kernel is only called 
   with range 0..1 for u as done in hydra or density !! 
   Call with mode 0 to calculate dwk and wk
   Call with mode -1 to calculate only wk
   Call with mode +1 to calculate only dwk */

inline void kernel_main(double u, double hinv3, double hinv4, 
        double *wk, double *dwk, int mode)
{
#ifdef STANDARD_KERNEL /* cubic spline */
  if(u < 0.5)
    {
      if(mode >= 0) 
          *dwk = u * (18.0 * u - 12.0);
      if(mode <= 0) 
          *wk = (1.0 + 6.0 * (u - 1.0) * u * u);
    }
  else
    {
      double t1 = (1.0 - u);
      double t2 = t1 * t1;
      if(mode >= 0) 
          *dwk = -6.0 * t2;
      if(mode <= 0) 
          *wk = 2.0 * t2 * t1;
    }
#endif /* STANDARD_KERNEL */

#ifdef QUINTIC_KERNEL
  double t1 = (1.0 - u);
  double t2 = t1 * t1;
  double t4 = t2 * t2;

  if(mode >= 0) 
      *dwk = -5.0 * t4;
  if(mode <= 0) 
      *wk = t4 * t1;

  if (u < 2.0/3.0)
    {
      t1 = (2.0/3.0 - u);
      t2 = t1 * t1;
      t4 = t2 * t2;
      if(mode >= 0) 
          *dwk += 30.0 * t4;
      if(mode <= 0) 
          *wk -= 6.0 * t4 * t1;
    }
  if (u < 1.0/3.0)
    {
      t1 = (1.0/3.0 - u);
      t2 = t1 * t1;
      t4 = t2 * t2;
      if(mode >= 0) 
          *dwk -= 75.0 * t4;
      if(mode <= 0) 
          *wk += 15.0 * t4 * t1;
    }
#endif /* QUINTIC_KERNEL */

#ifdef WENDLAND_C4_KERNEL /* Dehnen & Aly 2012 */
#if defined(THREEDIMS) || defined(TWODIMS)
  double t1 = (1-u), t5 = t1*t1*t1*t1*t1, u2 = u*u;

  if(mode >= 0) 
      *dwk = -280.0/3.0*t5 * u2  - 56.0/3.0*u*t5 ;  
  if(mode <= 0) 
      *wk = t1*t5*(1.0+6.0*u+35.0/3.0*u2);
#else /* ONEDIM */
  double t1 = (1-u), t4 = t1*t1*t1*t1, u2 = u*u;

  if(mode >= 0) 
      *dwk = -14.0 * u * t4 - 56.0 * u2 * t4;
  if(mode <= 0) 
      *wk =  t1*t4*(1.0+5.0*u+8.0*u2);
#endif /* ONEDIM */

#endif /* WENDLAND_C4_KERNEL */

#ifdef WENDLAND_C6_KERNEL /* Dehnen & Aly 2012 */
#if defined(THREEDIMS) || defined(TWODIMS)
  double t1 = (1-u), t7 = t1*t1*t1*t1*t1*t1*t1, u2 = u*u;

  if(mode >= 0) 
      *dwk = -22 * t7 * u * (16 * u2 + 7 * u + 1);
  if(mode <= 0) 
      *wk = t1 * t7 * (1 + 8 * u + 25 * u2 + 32 * u2 * u);
#else /* ONEDIM */
  double t1 = (1-u), t6 = t1*t1*t1*t1*t1*t1, u2 = u*u;

  if(mode >= 0) 
      *dwk = -6 * t6 * u * (35 * u2 + 18 * u + 3);
  if(mode <= 0) 
      *wk =  t6 * t1 * (1 + 7 * u + 19 * u2 + 21 * u2 * u);
#endif /* ONEDIM */

#endif /* WENDLAND_C6_KERNEL */


  if(mode >= 0) 
      *dwk *= NORM * hinv4;
  if(mode <= 0) 
      *wk *= NORM * hinv3;

  return;
}

#endif
