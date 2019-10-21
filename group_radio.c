#include "allvars.h"

#if defined(GROUPSPEC) || defined(OUTPUTGROUPLUM)
double group_luminosity( double nu, long index) {

    long p;
    double F;
    struct group_properties *g;
#ifdef GROUPLUMDENSITYWEIGHTED
    double d = 0;
#endif

    g = &Gprops[index];
    p = g->Head;
    F = 0;
    while( p >= 0 ) {

        if ( P[p].Type == 0 ) {

#ifdef GROUPLUMDENSITYWEIGHTED
            F += get_particle_radio( p , nu ) * SphP[p].Density;
            d += SphP[p].Density;
#else
            F += get_particle_radio( p, nu);
#endif
        }
        p = FoFNext[p];

    }
#ifdef GROUPLUMDENSITYWEIGHTED
    F /= d;
#endif

    return F;

}
#endif

#if defined(GROUPSPEC)
void group_flux( double nu, long index, double *flux, double *flux_nosr ) {

    double L;
    struct group_properties *g;
    double size;

    g = &Gprops[index];

    L = group_luminosity( nu, index );
    size = get_group_size( g );
    //*flux_nosr =  *flux = L;
    *flux_nosr = L / ( 4.0 * PI * SQR( LumDis * g2c.cm ) );
    *flux = *flux_nosr / ( SQR( size * 2  / ComDis ) );

}
#endif

void group_spectrum() {

#ifdef GROUPSPEC
    int vN, i, index;
    double v, *flux, *flux_nosr, vmin, vmax, dv;

    FILE *fd1, *fd2;

    put_header( "group spectrum" );

    vN = All.FreqN;
    vmin = All.FreqMin;
    vmax = All.FreqMax;

    mymalloc1( flux, sizeof(double) * vN );
    mymalloc1( flux_nosr, sizeof(double) * vN );

    dv = log( vmax/vmin ) / ( vN-1 );

    create_dir( "%sSpectrum", GroupDir );

    fd1 = myfopen("w", "%sSpectrum/Spec_%03i.dat",
            GroupDir, SnapIndex );

    fd2 = myfopen( "w", "%sSpectrum/Spec_%03i_nosr.dat",
            GroupDir, SnapIndex );

    mymalloc1( flux, sizeof(double) * vN );
    mymalloc1( flux_nosr, sizeof(double) * vN );

    fprintf( fd1, "0  0  " );
    fprintf( fd2, "0  0  " );

    for ( i=0; i<vN; i++ ) {
        v = exp(log(vmin) + i * dv);
        fprintf( fd1, "%g  ", v/1e6*Time );
        fprintf( fd2, "%g  ", v/1e6*Time );
    }
    fprintf( fd1, "\n" );
    fprintf( fd2, "\n" );

    for ( index=0; index<Ngroup; index++ ) {

        if ( !group_present( index ) )
            break;

        for ( i=0; i<vN; i++ ) {
            v = exp(log(vmin) + i * dv);
            group_flux( v, index, flux+i, flux_nosr+i );
        }

        fprintf( fd1, "%i  %g  ", index, Gprops[index].mass );
        fprintf( fd2, "%i  %g  ", index, Gprops[index].mass );

        for ( i=0; i<vN; i++ ) {
            fprintf( fd1, "%g  ", flux[i] );
            fprintf( fd2, "%g  ", flux_nosr[i] );
        }

        fprintf( fd1, "\n" );
        fprintf( fd2, "\n" );

    }

    fclose( fd1 );
    fclose( fd2 );

    myfree( flux );
    myfree( flux_nosr );
    put_end();

#endif
}

void group_rad_profile() {

#ifdef GROUPRADPROFILE
    int RN, gi, i, *n;
    double dlogr, rmin, rmax, *I, r, rr, fac_to_arcmin, fac_to_flux, area;
    group_properties g;
    long p;
    FILE *fd;

    put_header( "group rad profile" );
    RN = All.GroupRadProfileRN;
    mymalloc2( I, sizeof(double) * RN );
    mymalloc2( n, sizeof(int) * RN );
    //dlogr = log(rmax/rmin) / (RN-1);
    fd = myfopen( "w", "%s/GroupRadProfile_%03i.dat", GroupDir, SnapIndex );

    rr = 0.05;
    fac_to_flux = 1 / ( 4*PI * SQR(LumDis*g2c.cm) );
    fac_to_arcmin = SQR(1 / PI * 180 * 60);
    dlogr = log(1/rr) / (RN-1);

    for(i=0; i<RN; i++)
         fprintf( fd, "%g ", rr*exp(i*dlogr) );
    fprintf( fd, "\n" );

    for( gi=0; gi<Ngroup; gi++ ) {
        if ( !group_present(gi) )
            continue;

        g = Gprops[gi];
        p = g.Head;

        for( i=0; i<RN; i++ )
            I[i] = 0;
        rmax = g.vr200;;
        rmin = rr * rmax;
    
        while( p>=0 ) {
    
            if ( P[p].Type ) {
                p = FoFNext[p];
                continue;
            }

            for(i=0, r=0; i<3; i++)
                r += SQR( PERIODIC_HALF( P[p].Pos[i]-g.cm[i] ) );
            r = sqrt(r);

            if ( r > rmax ) {
                p = FoFNext[p];
                continue;
            }
    
            if ( r<rmin )
                i = 0;
            else
                i = log(r/rmin) / dlogr + 1;

            if ( i<=RN-1 ) {
                I[i] += get_particle_radio( p, All.GroupRadProfileFreq );
                n[i] ++;
            }
            p = FoFNext[p];
        }

        for( i=0; i<RN; i++ ) {
            if ( i==0 )
                area = PI * SQR(rmin);
            else
                area  = PI * (SQR(rmin*exp(i*dlogr))-SQR(rmin*exp((i-1)*dlogr)));
            area = area / SQR(ComDis) * fac_to_arcmin;
            I[i] = I[i] / area * fac_to_flux;
        }

        for(i=0; i<RN; i++)\
            fprintf( fd, "%g ", I[i] );
        fprintf( fd, "\n" );

        for(i=0; i<RN; i++)\
            fprintf( fd, "%i ", n[i] );
        fprintf( fd, "\n" );

    }

    fclose( fd );
    myfree( I );
    myfree( n );

    put_end();

#endif
}
#ifdef GROUPELECSPEC
double particle_f( SphParticleData *part, double p ) {

    double r;

    if ( p < part->CRE_qmin || p > part->CRE_qmax )
        return 0;

    r = part->CRE_C * part->Density / guc.m_e  * pow( p, -part->CRE_Alpha );
    r /= CUBE( g2c.cm );

    return r;

}
#endif

void group_electron_spectrum() {

#ifdef GROUPELECSPEC
    long p;
    int i, qn, index;
    double dlogq, qmin, qmax, *f, q, rho, qmax_max, qmin_min;


    FILE *fd;
#ifdef GROUP_ELECTRON_SPECTRUM_DEBUG
    FILE *fd2;
#endif

    qn = All.QNum;
    qmin = All.QMin;
    qmax = All.QMax;

    dlogq = log( qmax/qmin ) / ( qn-1 );

    create_dir(  "%sElectronSpectrum", GroupDir );
    fd = myfopen( "w", "%sElectronSpectrum/EleSpec_%03i.dat", GroupDir, SnapIndex );

    fprintf( fd, "0  0  " );
    for ( i=0; i<qn; i++ ) {
        q = exp(log(qmin) + i * dlogq);
        fprintf( fd, "%g  ", q );
    }
    fprintf( fd, "\n" );

    mymalloc2( f, sizeof(double) * qn );

    for ( index=0; index<Ngroup; index++ ) {

        if ( !group_present( index ) )
            break;

        p = Gprops[index].Head;
        rho = 0;
        qmax_max = 0;
        qmin_min = 1e10;

        for ( i=0; i<qn; i++ )
            f[i] = 0;

#ifdef GROUP_ELECTRON_SPECTRUM_DEBUG
        fd2 = myfopen( "w", "%sElectronSpectrum/PartInfo_%03i_%04i.dat", GroupDir, SnapIndex, index );

        fprintf( fd2, "%.2f %10.3e %10.3e 0 0 0\n",
                Redshift, SofteningTable[0] * g2c.cm,
                get_group_size( &Gprops[index] ) * g2c.cm );
#endif

        while( p >= 0 ) {
            if ( P[p].Type == 0 ) {
                for( i=0; i<qn; i++ ) {
                    q = log( qmin ) + i * dlogq;
                    q = exp(q);
#ifdef GROUPELECDENSITYWEIGHTED
                        f[i] += particle_f( &SphP[p], q ) * SphP[p].Density;
                        rho += SphP[p].Density;
#else
                        f[i] += particle_f( &SphP[p], q );
                        rho ++;
#endif
                }
                vmax2( qmax_max, SphP[p].CRE_qmax );
                vmin20( qmin_min, SphP[p].CRE_qmin );
#ifdef GROUP_ELECTRON_SPECTRUM_DEBUG
                fprintf( fd2, "%10.3e %6.2f %7.3f %8.0f %10.3e %10.3e\n",
                    SphP[p].CRE_C * SphP[p].Density / guc.m_e / CUBE( g2c.cm ),
                    SphP[p].CRE_Alpha,
                    SphP[p].CRE_qmin,
                    SphP[p].CRE_qmax,
                    SphP[p].Density,
                    get_B( p )
                    );
#endif
            }
            p = FoFNext[p];
        }
#ifdef GROUP_ELECTRON_SPECTRUM_DEBUG
        fclose( fd2 );
#endif

        //printf( "%li\n", Gprops[index].npart[0] );
        for ( i=0; i<qn; i++ ) {
            f[i] /= rho;
        }

        printf( "[%i] qmax_max: %g, qmin_min: %g\n", index, qmax_max, qmin_min );

        fprintf( fd, "%i  %g  ", index, Gprops[index].mass );
        for ( i=0; i<qn; i++ ) {
            fprintf( fd, "%g  ", f[i] );
        }
        fprintf( fd, "\n" );
    }

    myfree( f );

    fclose( fd );

#endif
}
