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
    /*
    double *lumgrid, *num, L;
    int j;
    */
    FILE *fd1, *fd2;

    put_header( "group spectrum" );

    vN = All.NuNum;
    vmin = All.NuMin * 1e6 /Time;
    vmax = All.NuMax * 1e6 /Time;

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
        fprintf( fd1, "%g  ", v );
        fprintf( fd2, "%g  ", v );
    }
    fprintf( fd1, "\n" );
    fprintf( fd2, "\n" );

    for ( index=0; index<Ngroups; index++ ) {

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

    for ( index=0; index<Ngroups; index++ ) {

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

/*
void group_spectrum_index() {

    double *spec, *v, vmin, vmax, dv, cov00, cov01, cov11, c0,
           *spec_index, *spec_index_err,L, dL, *mass;
    int vN, k, PicS, PicS2, ii, jj, i, j, x, y, xo, yo, flag,
        index, p;
    struct group_properties g;
    char buf[100],
         *spec_index_str="SpectrumIndex",
         *spec_index_err_str="SpectrumIndex_err";

    writelog( "group spectrum index... \n" );
    vN = All.NuNum;
    vmin = All.NuMin;
    vmax = All.NuMax;
    PicS = All.PicSize;
    PicS2 = SQR( PicS );
    x = proj_i;
    y = proj_j;
    xo = PicS / 2;
    yo = PicS / 2;

    dv = log10( vmax/vmin) / vN;

    create_dir( "%s%s", GroupDir, spec_index_str );

    mymalloc1( v, sizeof( double ) * vN );
    mymalloc1( spec, sizeof( double ) * vN * PicS2 );
    mymalloc1( spec_index, sizeof( double ) * PicS2 );
    mymalloc1( spec_index_err, sizeof( double ) * PicS2 );

    for ( i=0; i<vN; i++ )
        v[i] = log10(vmin) + i * dv;

    for ( index=0; index<Ngroups; index++ ) {

        if ( !group_present( index ) ) {
            break;
        }

        g = Gprops[index];
        memset( spec, 0, sizeof( double ) * vN * PicS2 );
        memset( spec_index, 0, sizeof( double ) * PicS2 );
        memset( spec_index_err, 0, sizeof( double ) * PicS2 );

        L = get_group_size( &g ) * 1.1;
        dL = 2 * L / PicS;

        p = g.Head;
        for ( i=0; i<g.Len; i++, p=FoFNext[p] ) {

            if ( P[p].Type != 0 )
                continue;

            ii = PERIODIC( P[p].Pos[x] - g.cm[x] ) / dL + xo;
            jj = PERIODIC( P[p].Pos[y] - g.cm[y] ) / dL + yo;

            check_picture_index( ii );
            check_picture_index( jj );

            for ( k=0; k<vN; k++ )
                spec[ii*PicS*vN + jj * vN + k] += get_particle_radio_index( p, k );
        }

        for ( i=0; i<vN * PicS2; i++ )
            if ( spec[i] != 0 )
                spec[i] = log10( spec[i] );

        for ( i=0; i<PicS; i++ )
            for ( j=0; j<PicS; j++ ) {
                flag = 0;
                for ( k=0; k<vN; k++ )
                    if ( spec[ i*PicS*vN + j*vN + k ] != 0 ){
                        flag = 1;
                        break;
                    }
                if ( flag ) {
                    gsl_fit_linear( v, 1, &spec[ i*PicS*vN + j*vN], 1,
                            vN, &c0, &spec_index[i*PicS + j],
                            &cov00, &cov01, &cov11,
                            &spec_index_err[i*PicS+j] );
                    if ( spec_index_err[i*PicS+j] > 1 ) {
                        spec_index_err[i*PicS+j] = 0;
                        spec_index[i*PicS+j] = 0;
                    }
                }
            }
        mass = g.mass_table;
        for ( i=0; i<6; i++ ) {
            mass[i] *= 1e10;
        }
        for ( i=0; i<6; i++ )
            img_props(i) = mass[i];
        img_xmin =  -L;
        img_xmax =  L;
        img_ymin =  -L;
        img_ymax =  L;

        image.img = spec_index;
        make_group_output_filename( buf, spec_index_str, index );
        write_img( buf );

        image.img = spec_index_err;
        sprintf( buf, "%s%s/%s_%03i_%04i_%c.dat",
            GroupDir, spec_index_str, spec_index_err_str,
            SnapIndex, index, Sproj );
        write_img( buf );

    }

    myfree( v );
    myfree( spec );
    myfree( spec_index );
    myfree( spec_index_err );
    reset_img();
    writelog( "group spectrum index... done.\n" );

}
*/
