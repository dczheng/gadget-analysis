#include "allvars.h"

double particle_df( double p, void *params ) {
    double *pa;
    pa = params;
    if ( p < pa[2] || p > pa[3] )
        return 0;
    return pa[0] * pow( p, -pa[1] );
}

double particle_radio( double nu,  SphParticleData *part ) {

    double r, params[4], B;

    params[0] = part->CRE_C;
    if ( part->CRE_qmin <= 0 )
        return 0;

    params[0] = params[0] * part->Density / guc.m_e;
    params[0] /= CUBE( g2c.cm );
    params[0] /= Time3;

    params[1] = part->CRE_Alpha;
    params[2] = part->CRE_qmin;
    params[3] = part->CRE_qmax;

    if (  params[ 0] *
            params[1] *
            params[2] *
            params[3] == 0 )
        return 0;

    B =  get_B( (part-SphP) );
    r = radio( &particle_df, params, B, nu, params[2], params[3], 1e-2 );
    return r;

}

void compute_particle_radio() {

#if defined(RAD) && !defined(RADLARMOR)
    double numin, numax, dlognu, nu, r_larmor, r;
    int nuN, j, *flags, *flags_local, sig;
    long i, dsig, idx, flags_sum, N, Ndo;

#ifdef PART_RADIO_TEST2
    double err=0, err_max=0, err_mean=0;
#endif

    writelog( "Start compute particle radio ... \n" )

    init_compute_F();

    mymalloc1_shared( PartRad, sizeof(double)*N_Gas*All.NuNum, sizeof(double), MpiWin_PartRad );

#ifdef TABF
    writelog( "use tab F\n" );
    init_tab_F();
#else
    writelog( "don't use tab F\n" );
#endif

    mymalloc2( flags, sizeof(int) * N_Gas );
    mymalloc2( flags_local, sizeof(int) * N_Gas );

    nuN = All.NuNum;
    numin = All.NuMin;
    numax = All.NuMax;
    sig = 10;
    dlognu = log( numax/numin ) / ( nuN - 1 );

    writelog( "numin: %g, numax: %g, sig: %i\n", numin, numax, sig );
    dsig  = N_Gas / sig;

    writelog( "first\n" );
    for( i=0; i<N_Gas; i++ ) {

        if ( i % dsig == 0 )
            writelog( "[%10li] [%10li] [%5.1f%%]\n",\
            i, N_Gas, ((double)(i)) / N_Gas * 100 );

        if ( i % NTask_Local != ThisTask_Local )
            continue;

        flags_local[i] = 1;
        //flags[i] = 1;
        for( j=0; j<nuN; j++ ) {
            nu = exp( log(numin) + j * dlognu );
            r_larmor = particle_radio_larmor( nu, &SphP[i] );
            if (r_larmor>0) {
                flags_local[i] =0;
                //flags[i] = 0;
                break;
            }
        }
    }

    MPI_Allreduce( flags_local, flags, N_Gas, MPI_INT, MPI_SUM, MpiComm_Local );
    myfree( flags_local );
    do_sync_local( "" );

    for( i=0, flags_sum=0; i<N_Gas; i++ ) {
        flags_sum += flags[i];
    }

    Ndo = N_Gas-flags_sum;
    //writelog( "second [%li] %li \n", N_Gas, flags_sum );
    writelog( "second [%li][%.2f%%]\n", Ndo, ((double)Ndo)/N_Gas * 100 );

    dsig = Ndo / sig;
    N = 0;
    for( i=0; i<N_Gas; i++ ) {

        if ( flags[i] )
            continue;

#ifdef PART_RADIO_TEST2
        if ( ThisTask_Local )
            continue;
#endif

        for( j=0; j<nuN; j++ ) {
            idx = i * nuN + j;

#ifndef PART_RADIO_TEST2
            if ( idx % NTask_Local != ThisTask_Local )
                continue;
#endif
            nu = exp( log(numin) + j * dlognu );
#if defined(PART_RADIO_TEST2) || defined(PART_RADIO_TEST3)
#ifdef PART_RADIO_TEST3
            if ( !ThisTask_Local )
#endif
                printf(
                        "--\n[%s], (%li, %i) M: %g, rho: %g, h: %g, B: %g, nu: %g,"
                        " c: %g, a: %g, qmin: %g, qmax: %g\n",
                    __FUNCTION__, i, j,
                    P[i].Mass, SphP[i].Density, SphP[i].Hsml,
                    get_B(i), nu,
                    SphP[i].CRE_C,
                    SphP[i].CRE_Alpha,
                    SphP[i].CRE_qmin,
                    SphP[i].CRE_qmax
                    );
#endif
            r = particle_radio( nu, &SphP[i] );
#ifdef PART_RADIO_TEST2
            r_larmor = particle_radio_larmor( nu, &SphP[i] );
            if ( r_larmor ) {
                err = (r_larmor-r)/r * 100;
                printf( "r: %g, r_larmor: %g, err: %.2f%%\n",\
                    r, r_larmor, err );
            }
            else {
                err = 0;
            }

            if ( err>err_max )
                err_max = err;
            err_mean += err;
#endif
            PartRad[idx] = r;
        } // for j

        if ( N % dsig == 0 )
            writelog( "[%10li] [%10li] [%5.1f%%]\n",\
                N, Ndo, ((double)(N)) / Ndo * 100 );
        N++;
    } // for i
#ifdef PART_RADIO_TEST2
    if ( !ThisTask_Local ) {
        err_mean /= Ndo;
        printf( "err_max: %.2f%%, err_mean: %.2f%%\n",
                    err_max, err_mean );
    }
#endif

    for( i=0; i<N_Gas; i++ ) {
        if ( !flags[i] )
            continue;
        for( j=0; j<nuN; j++ ) {
            idx = i*nuN + j;
            PartRad[idx] = 0; 
        }
    }

    myfree( flags );

#ifdef TABF
    free_tab_F();
#endif
    do_sync_local( "" );
#endif
}

#if defined(RADLARMOR) || defined(RAD)
double particle_radio_larmor( double nu, SphParticleData *part ) {

    double vL, B, a, pmin, pmax, C, t, r;
    long idx;
    idx = part - SphP;

    B = get_B( idx );

    if(B==0)  return 0;

    vL = B / ( 2*PI ) * cuc.e_mec;

    C = SphP[idx].CRE_C * SphP[idx].Density / (guc.m_e * CUBE(g2c.cm) * Time3);
    a = (SphP[idx].CRE_Alpha-1)/2.0;
    pmin = SphP[idx].CRE_qmin;
    pmax = SphP[idx].CRE_qmax;

    if ( nu<vL )
        return 0;
    t = sqrt(nu/vL-1);
    if ( t<pmin || t>pmax )
        r = 0;
    else
        r = 2.0/3.0 * cuc.sigma_t * C * cuc.c * SQR(B) / ( 8*PI * vL ) * pow( nu/vL-1, -a ); 
    return r;
}
#endif


double get_particle_radio( long p, double nu ) {

/*
nu: Mhz
*/
    double V, r;
    V = get_V(p);
#if defined(RAD) && !defined(RADLARMOR)
    double nu_x, dlognu;
    int nu_i;
    dlognu = log( All.NuMax/ All.NuMin  ) / ( All.NuNum -1 );
    nu_x = log( nu / All.NuMin ) / dlognu;
    nu_i = nu_x;
    nu_x -= nu_i;
    if ( nu_i >= All.NuNum-1  )
        r =  PartRad[p*All.NuNum+nu_i+All.NuNum-1] ;
    if ( nu_i < 0 )
        r =  PartRad[p*All.NuNum+nu_i] ;

    r =  PartRad[p*All.NuNum+nu_i] * ( 1-nu_x  ) + 
            PartRad[p*All.NuNum+nu_i+1] * nu_x;
#endif

#ifdef RADLARMOR
    r =  particle_radio_larmor( nu, &SphP[p] );
#endif

    return r*V;
}

void test_part_radio() {

#ifdef PART_RADIO_TEST

    HubbleParam = 0.68;
    Redshift = 0.104325;
    Time = (1/(Redshift+1));
    Omega0 = 0.302;
    OmegaLambda = 0.698;
    OmegaBaryon = 0.0471;
    set_units();
    compute_cosmo_quantities();

    mymalloc2( SphP, sizeof(SphParticleData) );
    mymalloc2( P, sizeof(ParticleData) );
    double nu;
    double r, r_larmor;

    init_compute_F();

#ifdef TABF
    init_tab_F();
#endif

    if ( !ThisTask_Local ) {
        /*
        P[0].Mass =  0.0785971;
        SphP[0].Hsml = 2.2638;
        SphP[0].Density = 0.052087;
        SphP[0].B[0] = SphP[0].B[1] = 0;
        SphP[0].B[2] = 2.7017e-5;
        SphP[0].CRE_C = 4.12802e-13; 
        SphP[0].CRE_Alpha = 2.42518;
        SphP[0].CRE_qmin = 221.434; 
        SphP[0].CRE_qmax = 1.55059e+06; 
        */

        P[0].Mass =  0.0785971;
        SphP[0].Density = 1.95961e-08;
        SphP[0].Hsml = 311.492;
        SphP[0].B[0] = SphP[0].B[1] = 0;
        SphP[0].B[2] = 9.33349e-10;
        SphP[0].CRE_C = 4.29783e-12; 
        SphP[0].CRE_Alpha = 2.92749;
        SphP[0].CRE_qmin = 0.19628; 
        SphP[0].CRE_qmax =96857.1; 

        nu = 3.40409e+08;

        r = particle_radio( nu,  &SphP[0]  );
        r_larmor = particle_radio_larmor( nu,  &SphP[0]  );
    
        printf( "r: %g, r_larmor: %g, err: %g%%\n", r, r_larmor, (r-r_larmor)/r*100 );
        printf( "f1: %g\n", sqrt(3)/4.0 * PI * CUBE(cuc.e) / cuc.mec2 );
        printf( "f2: %g\n", 3.0/16.0 * cuc.e_mec );
        printf( "vL: %g\n", 1e-6 / ( 2 * PI ) * cuc.e_mec);
        myfree( P );
        myfree( SphP );
    }

#ifdef TABF
    free_tab_F();
#endif
    do_sync( "" );
    if ( !ThisTask_Local )
        endruns("test_part_radio");
#endif

}


void free_particle_radio() {
#if defined(RAD) && !defined(RADLARMOR)
    myfree_shared( MpiWin_PartRad );
#endif
}
