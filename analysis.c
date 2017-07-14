#include "allvars.h"

void magnetic_field_analysis() {
     long i, index;
     double sum, b, bmax, min, max;
     sum = 0;
     index = 0;
     bmax = 0;
     for ( i=0; i<Particle[0].num; i++ ) {
         /*
         if( Particle[0].mag[i*3+0] != 0 ||
             Particle[0].mag[i*3+1] != 0 ||
             Particle[0].mag[i*3+2] != 0 )
         fprintf( stdout, "( %.2f, %.2f, %.2f ): %e %e %e\n",
                 Particle[0].pos[i*3+0],
                 Particle[0].pos[i*3+1],
                 Particle[0].pos[i*3+2],
                 Particle[0].mag[i*3+0],
                 Particle[0].mag[i*3+1],
                 Particle[0].mag[i*3+2] );
                 */
         b = pow( Particle[0].mag[i*3+0], 2 ) +
             pow( Particle[0].mag[i*3+1], 2 ) +
             pow( Particle[0].mag[i*3+2], 2 );
         if ( b>bmax ) {
             index =i;
             bmax = b;
         }
         sum += b;
     }
    fprintf( stdout, "max magnetic field: ( %.2f, %.2f, %.2f ): %e %e %e %e\n",
            Particle[0].pos[index*3+0],
            Particle[0].pos[index*3+1],
            Particle[0].pos[index*3+2],
            Particle[0].mag[index*3+0],
            Particle[0].mag[index*3+1],
            Particle[0].mag[index*3+2], sqrt( bmax ) );
     fprintf( stdout, "Total Magnetic Field: %e\n", sqrt( sum / Particle[0].num ) );
}

void density_analysis() {
     long i, index;
     double sum, min, max, rho;
     fputs( sep_str, stdout );
     fputs( "analyze density ...\n", stdout );
     sum = 0;
     index = 0;
     /*
     for ( i=0; i<Particle[0].num; i++ ) {
         rho = Particle[0].rho[i] * 1.989 * 1e43 / pow( 3.08e21,3 ); //* 150000 * 3.085e21;
         fprintf( stdout, "%e\n", rho );
         if ( i==10 ) break;
     }
     max = -1;
     min = 1000;
     for ( i=0; i<Particle[0].num; i++ ) {
        min = ( Particle[0].m[i]<min ) ? Particle[0].m[i] : min;
        max = ( Particle[0].m[i]>max ) ? Particle[0].m[i] : max;
     }
     fprintf( stdout, "min = %lf, max = %lf\n", min, max );
     min = 1e10;
     max = -1;
     for ( i=0; i<Particle[0].num; i++ ) {
        min = ( Particle[0].rho[i]<min ) ? Particle[0].rho[i] : min;
        max = ( Particle[0].rho[i]>max ) ? Particle[0].rho[i] : max;
     }
     fprintf( stdout, "min = %e, max = %e\n", min, max );
     */
    sum = 0;
     for ( i=0; i<Particle[0].num; i++ ) {
         sum += Particle[0].m[i];
     }
     fprintf( stdout, "m=%e\n", sum );
     fputs( sep_str, stdout );
}

int compare_for_sort_group_by_mass( const void *a, const void *b ) {
    return (( ( struct group_struct* )a )->Mass < ( ( struct group_struct* )b )->Mass ) ? 1 : -1;
}

void sort_group_by_mass() {
    qsort( (void*)group, TotNgroups, sizeof( struct group_struct ), &compare_for_sort_group_by_mass );
}

void show_group_info() {
    int i;
    for ( i=0; i<TotNgroups; i++ ) {
        if ( group[i].Mass >= 1e3 )
        fprintf( stdout, "mass: %10.3f, cm: ( %10.3f, %10.3f, %10.3f ), vel: ( %10.3f, %10.3f, %10.3f )\n",
                group[i].Mass, group[i].CM[0], group[i].CM[1], group[i].CM[2],
                group[i].Vel[0], group[i].Vel[1], group[i].Vel[2] );
    }
}

void group_analysis() {
    read_group();
    sort_group_by_mass();
    show_group_info();
    free_group();
}

void velocity_analysis() {
    int pt;
    long i;
    float vmax[3], v;
    pt = 0;
    vmax[0] = Particle[pt].vel[0];
    vmax[1] = Particle[pt].vel[1];
    vmax[2] = Particle[pt].vel[2];
    v = sqrt( pow( vmax[0], 2 ) + pow( vmax[1], 2 ) + pow( vmax[2], 2 ) );
    for ( i=0; i<Particle[pt].num; i++ ) {
        if ( sqrt( pow( Particle[pt].vel[ i*3+0 ], 2 ) +
                   pow( Particle[pt].vel[ i*3+1 ], 2 ) +
                   pow( Particle[pt].vel[ i*3+2 ], 2 ) ) > v ) {
            vmax[0] = Particle[pt].vel[ i*3+0 ];
            vmax[1] = Particle[pt].vel[ i*3+1 ];
            vmax[2] = Particle[pt].vel[ i*3+2 ];
            v = sqrt( pow( vmax[0], 2 ) +
                    pow( vmax[1], 2 ) +
                    pow( vmax[2], 2 ) );
        }
    }
    fprintf( stdout, "Vmax: %f ( %f, %f, %f )\n", v, vmax[0],
            vmax[1], vmax[2] );
}

double gamma_integrand( double t, void *params ) {
    double *x;
    x = ( double* ) params;
    return pow( t, *x-1 ) * exp( -t );
}

double gamma( double x ) {
    double epsabs, epsrel, abserr, result;
    size_t subinter;
    gsl_function F;
    epsabs = epsrel = 1e-8;
    subinter = 100000;
    gsl_integration_workspace *integration_workspace =
        gsl_integration_workspace_alloc( subinter );
    F.function = gamma_integrand;
    F.params = &x;
    gsl_integration_qagiu( &F, 0, epsabs, epsrel, subinter,
            integration_workspace, &result, &abserr );
    gsl_integration_workspace_free( integration_workspace );
    return result;
}

void analysis_radio() {
    double alpha_e, alpha, sigma_pp, sigma_t, me, mp, GeV, C, e;
    double Bcmb, alpha_v, XHe, Ecmb, Ce, Ae, Bc, C_phy, Eb, v, gamma_fac, q_phy;
    long i;
    alpha = 2.5;
    sigma_pp = 32 * ( 0.96 + exp( 4.4 - 2.4 * alpha ) ) * 1e-28 * 1e4; // cm^2
    me = 9.10938356e-29; // g
    mp = 1.672621898e-24; // g
    sigma_t = 6.6524587158e-29 * 1e4; // cm^2;
    C = 29979245800.0; // m/s
    GeV = 1.602176565e-19 * 1e9 * 1e7; // erg
    Bcmb = 3.24 * pow( 1+redshift, 2 ) * 1e-6; // G
    alpha_e = alpha + 1;
    alpha_v = alpha / 2.0;
    XHe = 0.24;
    e = 1.6021766208e-19; // C
    Ecmb = Bcmb * Bcmb / 8.0 / M_PI;
    v = 150000;
    gamma_fac = gamma( (3*alpha_e-1) / 12.0 ) *
            gamma( (3*alpha_e+7) / 12.0 ) *
            gamma( (alpha_e+5) / 4.0 ) /
            gamma( (alpha_e+7) / 4.0 );
    fprintf( stdout, "alpha = %f\nalpha_e = %f\nalpha_v = %e\n"
            "sigma_pp = %e\nsigma_t = %e\nC = %e\nGeV = %e\n"
            "me = %e\nmp = %e\nBcmb = %e\nXHe = %f\ne = %e\n"
            "Ecmb = %e\ngamma_fac = %e\n",
            alpha, alpha_e, alpha_v, sigma_pp, sigma_t, C, GeV, me, mp,
            Bcmb, XHe, e, Ecmb, gamma_fac );
    Bc = 2.0 * M_PI * pow( me, 3 ) * pow( C, 5 ) * v / ( 3 * e * pow( GeV, 2 ) );
    for ( i=0; i<Particle[0].num; i++ ) {
        if ( Particle[0].c0[i] == 0 ) {
            Particle[0].j[i] = 0;
            continue;
        }
        C_phy =  Particle[0].c0[i] * pow( Particle[0].rho[i], ( alpha - 1 ) * 0.33333 ) *
            Particle[0].rho[i] * 1.989e43 / pow( 3.085678e21, 3) / mp;
        //q_phy = Particle[0].q0[i] * pow( Particle[0].rho[i], 0.3333333 );
        Eb = ( pow( Particle[0].mag[i*3 + 0], 2 ) +
             pow( Particle[0].mag[ i*3 +1 ], 2 ) +
             pow( Particle[0].mag[ i*3 + 2 ], 2 ) ) / 8.0 / M_PI;
        Ce = pow( 16, 2-alpha_e ) / ( alpha_e-2 ) *
            sigma_pp * pow( me,2 ) * pow( C, 4 ) / ( sigma_t*GeV ) *
            C_phy * Particle[0].rho[i] * 1.989e43 / pow( 3.085678e21, 3) *
            Particle[0].elec[i] / ( 1-0.5*0.24 ) / ( Eb + Ecmb ) *
            pow( mp*C*C/GeV, alpha-1 );
        Ae = sqrt( 3.0 * M_PI ) / 32.0 / M_PI * Bc * pow( e, 3 ) / me / pow( C, 2 ) *
            ( alpha_e + 7.0/3.0 ) / ( alpha_e + 1 ) * gamma_fac;
        Particle[0].j[i] = 1e100 * Ae * Ce * pow( Eb/( pow(Bc,2) / 8.0 / M_PI ), ( alpha_v+1 ) / 2.0 );
        fprintf( stdout, "C_phy = %e, q_phy = %e, Eb = %e, Ce = %e, Bc = %e, Ae = %e, j = %e\n",
                C_phy, q_phy, Eb, Ce, Bc, Ae,  Particle[0].j[i] );
    }
    //plot_slice( 0, IO_J );
}
