#include "allvars.h"

#define make_group_output_filename( buf, nstr, group_index ) \
    sprintf( buf, "%s%s/%s_%03i_%04i_%c.dat",\
            GroupDir, nstr, nstr, SnapIndex, group_index, Sproj );

int group_present( long index ) {

    if ( Gprops[index].mass >= All.GroupMassMin) 
            return 1;
    else {

        return 0;
    }

}

/*

void group_particle_spectrum() {

    double numin, numax, dlognu, *nu;
    int nuN, index, i, len, l;
    long p, t;

    nuN = All.NuNum;
    t = nuN / 10;
    numin = All.NuMin;
    numax = All.NuMax;

    dlognu = log( numax/numin ) / ( nuN-1 );
    mymalloc1( nu, sizeof(double)*nuN );

    writelog( "compute group particle spectrum...\n" );

    put_sep0;

    for( i=0; i<nuN; i++ )
        nu[i] = exp( log(numin) + i * dlognu ) * 1e6;

    for ( index=0; index<Ngroups; index++ ) {

        if ( !group_present( index ) )
            break;

        len = Gprops[index].Len;
        t = len / 10;
        p = Gprops[index].Head;
        writelog( "group: %i ...\n", index );
        for ( l=0; l<len; l++, p = FoFNext[p]) {
            if ( l % t == 0 )
                writelog( "[%8i] [%8i] [%5.1f%%]\n", l, len, ((double)(l))/len * 100);
            if ( P[p].Type != 0 )
                continue;

            SphP[p].P = malloc( sizeof( double ) * nuN );
            for( i=0; i<nuN; i++ )
                SphP[p].P[i] = particle_radio( nu[i], p );

        }

    }

    myfree( nu );
    writelog( "compute group particle spectrum... done.\n" );
    put_sep0;

}

void free_group_particle_spectrum() {

    int index;
    long p;
    for ( index=0; index<Ngroups; index++ ) {

        if ( !group_present( index ) )
            continue;
        p = Gprops[index].Head;
        while( p >= 0 ) {

            if ( P[p].Type == 0 ) {
                free( SphP[p].P );
            }

            p = FoFNext[p];
        }

    }

}
*/

/*
void group_luminosity( int nu_index, long index, double *lumgrid, double *num ) {

    long p;
    int x[3], i, c, ii;
    struct group_properties *g;
    double L, dL, dL3;

    memset( lumgrid, 0, sizeof(double) * CUBE(All.GroupLumGrid) );
    memset( num, 0, sizeof(double) * CUBE(All.GroupLumGrid) );
    g = &Gprops[index];
    p = g->Head;
    c = All.GroupLumGrid / 2;

    L = 0;
    for(i=0; i<3; i++)
       L = vmax( L, g->size[i] );

    L *= 1.01;
    dL = 2 * L / All.GroupLumGrid;
    dL3 = CUBE(dL * Time * g2c.cm);

    while( p >= 0 ) {

        if ( P[p].Type != 0 ) {
            p = FoFNext[p];
            continue;
        }

        for( i=0; i<3; i++ ) {
            x[i] = PERIODIC_HALF( P[p].Pos[i] - g->cm[i] ) / dL + c;
            if ( x[i]<0 || x[i] > All.GroupLumGrid-1 ) {
                printf( "%s:warning!!!\n", __FUNCTION__ );
                p = FoFNext[p];
                continue;
            }
        }

        ii = get_index(x[0],x[1],x[2], All.GroupLumGrid);
        lumgrid[ii] += PartRad[ p * All.NuNum + nu_index] * SphP[p].Density
                        * dL3;
        num[ii] += SphP[p].Density; 
        p = FoFNext[p];
    }

    for(i=0; i<CUBE(All.GroupLumGrid); i++)
        if ( num[i] )
            lumgrid[i] /= num[i];

//#define GROUP_LUMINOSITY_TEST
#ifdef GROUP_LUMINOSITY_TEST
    FILE *fd;
    int ti, tj, tk;
    double sss;
    fd = fopen( "glum-test.dat", "w" );
    for( ti=0; ti<All.GroupLumGrid; ti++ ) {
        for( tj=0; tj<All.GroupLumGrid; tj++ ) {
            sss = 0;
            for( tk=0; tk<All.GroupLumGrid; tk++ )
                sss += lumgrid[ get_index(ti, tj, tk, All.GroupLumGrid) ];
            fprintf( fd, "%g ", sss );
        }
        fprintf( fd, "\n" );
    }
    fclose( fd );
    endruns( "group_luminosity-test" );
#endif

}

*/
double group_luminosity( int nu_index, long index ) {

    long p;
    double F;
    struct group_properties *g;

    g = &Gprops[index];
    p = g->Head;
    F = 0;

    while( p >= 0 ) {

        if ( P[p].Type == 0 )
            F += PartRad[ p * All.NuNum + nu_index];
        p = FoFNext[p];

    }

    return F;

}

void group_flux( int nu_index, long index, double *flux, double *flux_nosr ) {

    double L;
    struct group_properties *g;
    double size;

    g = &Gprops[index];

    L = group_luminosity( nu_index, index );
    size = get_group_size( g );
    *flux_nosr = L / ( 4.0 * PI * SQR( LumDis * g2c.cm ) );
    *flux = *flux_nosr / ( SQR( size * 2  / ComDis ) );

}



inline double get_group_size( struct group_properties *g ) {
    return vmax( g->size[proj_i], g->size[proj_j] );
}

void group_spectrum() {

    int vN, i, index;
    double v, *flux, *flux_nosr, vmin, vmax, dv;
    /*
    double *lumgrid, *num, L;
    int j;
    */
    FILE *fd1, *fd2;

    writelog( "group spectrum ...\n" );

    vN = All.NuNum;
    vmin = All.NuMin;
    vmax = All.NuMax;

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

/*
    mymalloc1( lumgrid, sizeof(double) * CUBE(All.GroupLumGrid) );
    mymalloc1( num, sizeof(double) * CUBE(All.GroupLumGrid) );
    */

    for ( index=0; index<Ngroups; index++ ) {

        if ( !group_present( index ) )
            break;


/*
        L = get_group_size( &Gprops[index] );
        memset( flux_nosr, 0, sizeof(double) *vN );
*/

        for ( i=0; i<vN; i++ ) {

/*
            group_luminosity( i, index, lumgrid, num );

            for(j=0; j<CUBE(All.GroupLumGrid); j++) {
                flux_nosr[i] += lumgrid[j];
            }

            flux_nosr[i] /= ( 4.0 * PI * SQR( LumDis * g2c.cm ) );
            flux[i] = flux_nosr[i] / ( SQR( L * 2  / ComDis ) );
*/
            group_flux( i, index, flux+i, flux_nosr+i );
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
    /*
    myfree( lumgrid );
    myfree( num );
    */

    writelog( "group spectrum ... done.\n" );

}


double particle_f( SphParticleData *part, double p ) {

    double r;

    if ( p < part->CRE_qmin || p > part->CRE_qmax )
        return 0;

    r = part->CRE_C * part->Density / guc.m_e  * pow( p, -part->CRE_Alpha );
    r /= CUBE( g2c.cm );

    return r;

}

void group_electron_spectrum() {

    long p;
    int i, qn, index;
    double dlogq, qmin, qmax, *f, q, rho, qmax_max, qmin_min;
    FILE *fd, *fd2;

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

        fd2 = myfopen( "w", "%sElectronSpectrum/PartInfo_%03i_%04i.dat", GroupDir, SnapIndex, index );

        fprintf( fd2, "%.2f %10.3e %10.3e 0 0 0\n",
                Redshift, SofteningTable[0] * g2c.cm,
                get_group_size( &Gprops[index] ) * g2c.cm );

        while( p >= 0 ) {
            if ( P[p].Type == 0 ) {
                for( i=0; i<qn; i++ ) {
                    q = log( qmin ) + i * dlogq;
                    q = exp(q);
                    //f[i] += particle_f( &SphP[p], q ) * SphP[p].Density;
                    f[i] += particle_f( &SphP[p], q );
                    rho += SphP[p].Density;
                }
                vmax2( qmax_max, SphP[p].CRE_qmax );
                vmin20( qmin_min, SphP[p].CRE_qmin );
                fprintf( fd2, "%10.3e %6.2f %7.3f %8.0f %10.3e %10.3e\n",
                    SphP[p].CRE_C * SphP[p].Density / guc.m_e / CUBE( g2c.cm ),
                    SphP[p].CRE_Alpha,
                    SphP[p].CRE_qmin,
                    SphP[p].CRE_qmax,
                    SphP[p].Density,
                    get_B( p )
                    );
            }
            p = FoFNext[p];
        }
        fclose( fd2 );

        //printf( "%li\n", Gprops[index].npart[0] );
        for ( i=0; i<qn; i++ ) {
            //f[i] /= rho;
            f[i] /= Gprops[index].npart[0];
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

}

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
                spec[ii*PicS*vN + jj * vN + k] += PartRad[ p * All.NuNum + k ];
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
        make_group_output_filename( buf, spec_index_str, index )
        write_img( buf, spec_index_str );

        image.img = spec_index_err;
        sprintf( buf, "%s%s/%s_%03i_%04i_%c.dat",
            GroupDir, spec_index_str, spec_index_err_str,
            SnapIndex, index, Sproj );
        write_img( buf, spec_index_err_str );

    }

    myfree( v );
    myfree( spec );
    myfree( spec_index );
    myfree( spec_index_err );
    reset_img();
    writelog( "group spectrum index... done.\n" );
    put_sep0;

}

int group_filed_present( enum group_fields blk ) {

    switch( blk ) {
        case GROUP_DENS:
                return 1;
        case GROUP_TEMP:
            if ( All.GroupTemp )
                return 1;
            return 0;
        case GROUP_SFR:
            if ( All.GroupSfr )
                return 1;
            return 0;
        case GROUP_MAG:
            if ( All.GroupB )
                return 1;
            return 0;
        case GROUP_MACH:
            if ( All.GroupMach )
                return 1;
            return 0;
        case GROUP_CREN:
        case GROUP_CREC:
        case GROUP_CREE:
        case GROUP_CREALPHA:
        case GROUP_CREQMIN:
        case GROUP_CREQMAX:
            if ( All.GroupCre )
                return 1;
            return 0;
        case GROUP_RAD:
            if ( All.GroupRad )
                return 1;
            return 0;
        default:
            return 0;
    }
}

void get_group_filed_name( enum group_fields blk, char *buf ) {
    switch( blk ) {
        case GROUP_DENS:
            strcpy( buf, "Density" );
            break;
        case GROUP_TEMP:
            strcpy( buf, "Temperature" );
            break;
        case GROUP_SFR:
            strcpy( buf, "SFR" );
            break;
        case GROUP_MAG:
            strcpy( buf, "MagneticField" );
            break;
        case GROUP_MACH:
            strcpy( buf, "Mach" );
            break;
        case GROUP_CREN:
            strcpy( buf, "Cre_n" );
            break;
        case GROUP_CREC:
            strcpy( buf, "Cre_C" );
            break;
        case GROUP_CREE:
            strcpy( buf, "Cre_e" );
            break;
        case GROUP_CREALPHA:
            strcpy( buf, "Cre_Alpha" );
            break;
        case GROUP_CREQMIN:
            strcpy( buf, "Cre_qmin" );
            break;
        case GROUP_CREQMAX:
            strcpy( buf, "Cre_qmax" );
            break;
        case GROUP_RAD:
            strcpy( buf, "Radio" );
            break;
    }
}

void check_group_flag() {

    int i, flag;
    flag = 0;

    for( i=0; i<GROUP_FILED_NBLOCKS; i++ ) {

        if ( !group_filed_present(i) )
            continue;

        switch( i ) {
            case GROUP_SFR:
                if ( All.ReadSfr == 0 ) {
                    printf( "`ReadSfr` is required by `GroupSfr`.\n" );
                    flag = 1;
                }
                break;
            case GROUP_MAG:
                if ( All.ReadB == 0 ) {
                    printf( "`ReadB` is required by `GroupB`.\n" );
                    flag = 1;
                }
                break;
            case GROUP_MACH:
                if ( All.ReadMach == 0 ) {
                    printf( "`ReadMach` is required by `GroupMach`.\n" );
                    flag = 1;
                }
                break;
            case GROUP_CREC:
            case GROUP_CREN:
            case GROUP_CREE:
            case GROUP_CREQMIN:
            case GROUP_CREQMAX:
            case GROUP_CREALPHA:
                if ( All.ReadCre == 0 ) {
                    printf( "`ReadCre` is required by `GroupCren`.\n" );
                    flag = 1;
                }
                break;
            case GROUP_RAD:
                if ( All.ReadCre == 0  || All.ReadB == 0) {
                    printf( "`ReadCre` and `RaadB` is required by `GroupRad`.\n" );
                    flag = 1;
                }
                break;
        }
    }

    if ( flag )
        endrun( 20180926 );

}

void group_temp_profile() {

    long p;
    struct group_properties g;
    int g_index, i, RN, *num, ii, ngbnum, x, y, z, k;
    double  *T, *data, Rmin, Rmax, dlogR, r;
    char buf[100], buf1[100], *dn="TempProfile";
    FILE *fd;

    RN = All.GroupTempProfileRN;
    Rmin = All.GroupTempProfileRmin;
    Rmax = All.GroupTempProfileRmax;
    dlogR = log10( Rmax/Rmin ) / ( RN-1 );

    mymalloc2( T,  sizeof(double) * RN );
    mymalloc2( num,  sizeof(int) * RN );
    mymalloc1( data, sizeof(double) * N_Gas * 3 );
    mymalloc1( Ngblist, NumPart * sizeof(long) );

    writelog( "Temperature profile ....\n" );
    create_dir( "%s%s", GroupDir, dn );

    reset_img();
    for ( g_index=0; g_index<Ngroups; g_index++ ) {

        if ( !group_present( g_index ) )
            break;

        g = Gprops[g_index];
        ngbnum = ngb( g.cm, Rmax*0.9 );
        /*
        ngbnum = 0;
        for ( p=0; p<N_Gas; p++ ) {
            for ( k=0, r=0; k<3; k++ )
                r += SQR( PERIODIC_HALF( P[p].Pos[k] - g.cm[k] ) );
            r = sqrt(r);
            if ( r < Rmax*0.9 )
                Ngblist[ngbnum++] = p;
        }
        */

        for( k=0; k<3; k++ ) {
            x = k;
            y = (k+1) % 3;
            z = (k+2) % 3;
            for( i=0; i<ngbnum; i++ ) {
                p = Ngblist[i];
                if ( P[p].Type != 0 )
                    continue;
                data[i*3] = PERIODIC_HALF( P[p].Pos[x] - g.cm[x] ) + Rmax;
                data[i*3+1] = PERIODIC_HALF( P[p].Pos[y] - g.cm[y] ) + Rmax;
                data[i*3+2] = SphP[p].Density / Time3 / RhoBaryon;
                //data[i*3+2] = SphP[p].Temp;
            }
            data_to_grid2d( data, image.img, ngbnum, All.PicSize,  2*Rmax );
            sprintf( buf, "%s%s/Density_%04i_%c.dat", GroupDir, dn, g_index, 'x'+z );
            sprintf( buf1, "Density_%04i_%c.dat", g_index, 'x'+k );
            write_img( buf, buf1 );
        }

        for( i=0; i<ngbnum; i++ ) {
            p = Ngblist[i];
            if ( P[p].Type != 0 )
                continue;
            for ( k=0, r=0; k<3; k++ )
                r += SQR( PERIODIC_HALF( P[p].Pos[k] - g.cm[k] ) );
            r = sqrt(r);
            ii = log10( r / Rmin ) / dlogR;
            if ( ii < 0 || ii > RN-1 )
                continue;
            //T[ii] += SphP[p].Temp;
            T[ii] += SphP[p].Density / Time3 / RhoBaryon;
            num[ii] ++;
        }
        for( i=0; i<RN; i++ )
            if ( num[i] > 0 )
                T[i] /= num[i];
        fd = myfopen( "w", "%s%s/TempProfile_%04i.dat", GroupDir, dn, g_index );
        for( i=0; i<RN; i++ )
            fprintf( fd, "%g %g %i\n", Rmin*pow(10, i*dlogR), T[i], num[i]  );
        fclose( fd );

    }

    myfree( T );
    myfree( num );
    myfree( data );
    myfree( Ngblist );

}

void group_temp_stack() {

    double Rmin, Rmax, dlogR, r;
    int RN, g_index, i, ii, mi, k, ngbnum;
    long p;
    struct group_properties g;
    FILE *fd;
#define MS 3
    double m[MS], *T[MS], *Terr[MS];
    int *num[MS]; 
    m[0] = 1e3;
    m[1] = 5e3;
    m[2] = 1e4;

    mymalloc1( Ngblist, NumPart * sizeof(long) );
    
    RN = All.GroupTempStackRN;
    Rmin = All.GroupTempStackRmin;
    Rmax = All.GroupTempStackRmax;
    dlogR = log10( Rmax/Rmin ) / (RN-1);

    for( i=0; i<MS; i++ ) {
        T[i] = malloc( sizeof(double) * RN );
        memset( T[i], 0, sizeof(double)*RN );
        Terr[i] = malloc( sizeof(double) * RN );
        memset( Terr[i], 0, sizeof(double)*RN );
        num[i] = malloc( sizeof(int) * RN );
        memset( num[i], 0, sizeof(int)*RN );
    }

    for ( g_index=0; g_index<Ngroups; g_index++ ) {
        g = Gprops[g_index];
        if ( g.mass < m[0] )
            continue;

        mi = 0;
        while( mi < MS-1 && g.mass > m[mi+1] ) mi ++;

        ngbnum = ngb( g.cm, All. GroupTempStackRmax * 1.1 );
        for( i=0; i<ngbnum; i++ ) {
            p = Ngblist[i];
            if ( P[p].Type != 0 )
                continue;

            for ( k=0, r=0; k<3; k++ )
                r += SQR( PERIODIC_HALF( P[p].Pos[k] - g.cm[k] ) );
            r = sqrt(r);
            ii = log10( r / Rmin ) /  dlogR;
            if ( ii < 0 || ii > RN-1 )
                continue;
            T[mi][ii] += SphP[p].Temp;
            num[mi][ii] ++;

        }

    /*
        for ( i=0,p=g.Head; i<g.Len; i++, p=FoFNext[p] ) {
            if ( P[p].Type != 0 )
                continue;

            for ( k=0, r=0; k<3; k++ )
                r += SQR( PERIODIC_HALF( P[p].Pos[k] - g.cm[k] ) );
            r = sqrt(r);
            ii = log10( r / Rmin ) /  dlogR;
            if ( ii < 0 || ii > RN-1 )
                continue;
            T[mi][ii] += SphP[p].Temp;
            num[mi][ii] ++;

        }
    */
    }

    for( mi=0; mi<MS; mi++ )
        for( i=0; i<RN; i++ )
            if ( num[mi][i] > 0 )
                T[mi][i] /= num[mi][i];

    for ( g_index=0; g_index<Ngroups; g_index++ ) {
        g = Gprops[g_index];
        if ( g.mass < m[0] )
            continue;
        mi = 0;
        while( mi < MS-1 && g.mass > m[mi+1] ) mi ++;

        ngbnum = ngb( g.cm, All. GroupTempStackRmax );
        for( i=0; i<ngbnum; i++ ) {
            p = Ngblist[i];
            if ( P[p].Type != 0 )
                continue;

            for ( k=0, r=0; k<3; k++ )
                r += SQR( PERIODIC_HALF( P[p].Pos[k] - g.cm[k] ) );
            r = sqrt(r);
            ii = log10( r / Rmin ) /  dlogR;
            if ( ii < 0 || ii > RN-1 )
                continue;
            Terr[mi][ii] += SQR( SphP[p].Temp - T[mi][ii] );

        }

/*
        for ( i=0,p=g.Head; i<g.Len; i++, p=FoFNext[p] ) {
            if ( P[p].Type != 0 )
                continue;

            for ( k=0, r=0; k<3; k++ )
                r += SQR( PERIODIC_HALF( P[p].Pos[k] - g.cm[k] ) );
            r = sqrt(r);
            ii = log10( r / Rmin ) /  dlogR;
            if ( ii < 0 || ii > RN-1 )
                continue;
            Terr[mi][ii] += SQR( SphP[p].Temp - T[mi][ii] );

        }
    */
    }

    for( mi=0; mi<MS; mi++ )
        for ( i=0; i<RN; i++ )
            if ( num[mi][i] > 1 )
                Terr[mi][i]  = sqrt(Terr[mi][i] / (num[mi][i]-1));


    create_dir( "%sGroupTempStack", OutputDir );
    fd = myfopen( "w", "%s/GroupTempStack/GroupTempStack_%03i.dat", OutputDir, SnapIndex );

    fprintf( fd, "%g ", Redshift );
    for( mi=0; mi<MS; mi++ )
        fprintf( fd, "%g 0 0 ", m[mi] );

    for( i=0; i<RN; i++ ) {
        fprintf( fd, "\n%g ", Rmin*pow( 10, i*dlogR ) );
        for( mi=0; mi<MS; mi++ )
            fprintf( fd, "%g %g %i ",
                T[mi][i], Terr[mi][i], num[mi][i] );
    }

    fclose( fd );


    for( mi=0; mi<MS; mi++ ) {
        free( T[mi] );
        free( Terr[mi] );
        free( num[mi] );
    }

    myfree( Ngblist );

}

void group_gas_ratio() {

    long p;
    struct group_properties g;
    int g_index, i;
    FILE *fd;

    double r_gas, r_condensed, r_diffuse, r_hot, r_warmhot, m_tot, m_gas;

    create_dir( "%sGroupGasRatio", OutputDir );
    fd = myfopen( "w", "%s/GroupGasRatio/GroupGasRatio_%03i.dat", OutputDir, SnapIndex );

    for ( g_index=0; g_index<Ngroups; g_index++ ) {
        if ( !group_present( g_index ) )
            break;
        g = Gprops[g_index];

        for( i=0,m_tot=0; i<6; i++ )
            m_tot += g.mass_table[i];
        m_gas = g.mass_table[0];
        r_gas = m_gas / m_tot;

        r_condensed = r_diffuse = r_hot = r_warmhot = 0; 
        p = g.Head;
        while( p >= 0 ) {
            if ( P[p].Type != 0 ) {
                p = FoFNext[p];
                continue;
            }

            if ( SphP[p].Temp >= 1e7 )
                r_hot += P[p].Mass;

            if ( SphP[p].Temp < 1e7 && SphP[p].Temp >= 1e5 )
                r_warmhot += P[p].Mass;

            if ( SphP[p].Temp < 1e5 && SphP[p].Density / Time3 / RhoBaryon < 1e3 )
                r_diffuse += P[p].Mass;

            if ( SphP[p].Temp < 1e5 && SphP[p].Density / Time3 / RhoBaryon >= 1e3 )
                r_condensed += P[p].Mass;

            p = FoFNext[p];
        }

        r_hot /= m_gas;
        r_warmhot /= m_gas; 
        r_condensed /= m_gas;
        r_diffuse /= m_gas;

        fprintf( fd, "%i %g %g %g %g %g\n", g_index,
            r_gas, r_hot, r_warmhot, r_condensed, r_diffuse );

    }

    fclose( fd );
}

void group_plot() {

    long p;
    struct group_properties g;
    double L, dL, *mass, B, PP, n, r200, dlognu, nu_x; 
    int *num, g_index, i, j, x, y,
         xo, yo, pic_index, ii, jj, nu_i;

    char buf[100], buf1[100];
    double *data[GROUP_FILED_NBLOCKS];

    x = proj_i;
    y = proj_j;
    xo = PicSize / 2;
    yo = PicSize / 2;

    dlognu = log( All.NuMax/ All.NuMin ) / ( All.NuNum );

    mymalloc2( num, PicSize2 * sizeof( int ) );

    for( i=0; i<GROUP_FILED_NBLOCKS; i++ ) {

        if ( !group_filed_present( i ) )
            continue;

        mymalloc1( data[i], PicSize2 * sizeof( double ) );
        get_group_filed_name( i, buf1 );
        create_dir( "%s%s", GroupDir, buf1 );

    }

    if ( All.GroupFixedSize ){
        if ( All.GroupSize > 0 ) {
            L = All.GroupSize;
        }
        else{
    
            L = 0;
            for ( g_index=0; g_index<Ngroups; g_index++ ) {
    
                if ( !group_present( g_index ) )
                    break;
                g = Gprops[g_index];
    
                for ( i=0,p=g.Head; i<g.Len; i++, p=FoFNext[p] ) {
                    if ( P[p].Type != 0 )
                        continue;
                vmax2( L, PERIODIC_HALF( P[p].Pos[x] - g.cm[x] ) * 1.01 );
                vmax2( L, PERIODIC_HALF( P[p].Pos[y] - g.cm[y] ) * 1.01 );
                }
            }
        }
    }

    for ( g_index=0; g_index<Ngroups; g_index++ ) {

        if ( !group_present( g_index ) )
            break;

        memset( num, 0, PicSize2 * sizeof( int ) );
        for( i=0; i<GROUP_FILED_NBLOCKS; i++ ) {

            if ( !group_filed_present(i) )
                continue;
            memset( data[i], 0, PicSize2 * sizeof( double ) );

        }

        writelog( "group plot[%c]: %i ...\n", Sproj, g_index );

//#define PRINT_GROUP_INFO
        g = Gprops[g_index];

#ifdef PRINT_GROUP_INFO
        writelog( "center of mass: %g %g %g\n",
                g.cm[0], g.cm[1], g.cm[2] );
#endif

        p = g.Head;
        mass = g.mass_table;
        r200 = g.vr200;

#ifdef PRINT_GROUP_INFO
        writelog( "npart: " );
        for ( i=0; i<6; i++ )
            writelog( "%li ", g.npart[i] );
        writelog( "\n" );

        writelog( "mass: " );
        for ( i=0; i<6; i++ ) {
            mass[i] *= 1e10;
            writelog( "%g ", mass[i] );
        }
#endif

        if ( !All.GroupFixedSize ) {
            L = 0;
            for ( i=0,p=g.Head; i<g.Len; i++, p=FoFNext[p] ) {
                if ( P[p].Type != 0 )
                    continue;
            vmax2( L, PERIODIC_HALF( P[p].Pos[x] - g.cm[x] ) * 1.1 );
            vmax2( L, PERIODIC_HALF( P[p].Pos[y] - g.cm[y] ) * 1.1 );
            }
        }

        //L = get_group_size( &g ) * 1.1;
        dL = 2 * L / PicSize;
#ifdef PRINT_GROUP_INFO
        writelog( "\n" );
        writelog( "L: %g, dL:%g\n", 2*L, dL );
#endif

        p = g.Head;
        for ( i=0; i<g.Len; i++, p=FoFNext[p] ) {
            if ( P[p].Type != 0 )
                continue;

            //printf( "[%li] Pos: ( %g, %g, %g )\n",
            //        p, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2] );
            //continue;

            ii = PERIODIC_HALF( P[p].Pos[x] - g.cm[x] ) / dL + xo;
            jj = PERIODIC_HALF( P[p].Pos[y] - g.cm[y] ) / dL + yo;

            //if ( ii < 0 || ii >= PicSize ||
            //        jj < 0 || jj >= PicSize )
            //        printf( "%i, %i\n", ii, jj );
            //
            //continue;

            check_picture_index( ii );
            check_picture_index( jj );
            pic_index= ii*PicSize + jj;

            num[ pic_index ] ++;

            if ( group_filed_present( GROUP_DENS ) )
                data[GROUP_DENS][pic_index] += SphP[p].Density;

            if ( group_filed_present( GROUP_TEMP ) ) {
                //printf( "%g\n", SphP[p].Temp );
                data[GROUP_TEMP][pic_index] += SphP[p].Temp * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_SFR ) )
                data[GROUP_SFR][pic_index] += SphP[p].sfr * SphP[p].Density;

            if ( group_filed_present( GROUP_MAG ) ) {
                B = sqrt( SQR( SphP[p].B[0] ) +
                          SQR( SphP[p].B[1] ) +
                          SQR( SphP[p].B[2] ) );
                B *= 1e6;  // convert G to muG
               // printf( "%g\n", B );
                data[GROUP_MAG][pic_index] += B * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_MACH ) )
                data[GROUP_MACH][pic_index] += SphP[p].MachNumber * SphP[p].Density;

            if ( group_filed_present( GROUP_CREN ) ) {
                n = SphP[p].CRE_n *  SphP[p].Density / guc.m_e;
                n /= CUBE( g2c.cm );
                data[GROUP_CREN][pic_index] += n * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_CREALPHA ) ) {
                data[GROUP_CREALPHA][pic_index] += SphP[p].CRE_Alpha * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_CREC ) ) {
                data[GROUP_CREC][pic_index] +=  SphP[p].CRE_C * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_CREE ) ) {
                /*
                e = SphP[p].CRE_e * SphP[p].Density;
                e /= g2c.erg / CUBE( g2c.cm );
                data[GROUP_CREE][pic_index] += e * SphP[p].Density;
                */
                data[GROUP_CREE][pic_index] += SphP[p].CRE_e / SphP[p].u * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_CREQMIN ) ) {
                data[GROUP_CREQMIN][pic_index] += SphP[p].CRE_qmin * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_CREQMAX ) ) {
                data[GROUP_CREQMAX][pic_index] += SphP[p].CRE_qmax * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_RAD ) ) {
                //PP = particle_radio( nu, p );
                nu_x = log(All.GroupRadFreq / All.NuMin) / dlognu;
                nu_i = nu_x;
                nu_x -= nu_i;
                if ( nu_i >= All.NuNum-1 )
                    PP = PartRad[p * All.NuNum + All.NuNum-1];
                else if ( nu_i <= 0 )
                    PP = PartRad[p * All.NuNum];
                else  {
                    PP = PartRad[p * All.NuNum+nu_i] * ( 1-nu_x ) +
                         PartRad[p * All.NuNum+nu_i+1] * nu_x;
                }
                data[GROUP_RAD][pic_index] += PP * SphP[p].Density;
            }

        }

        for ( i=0; i<PicSize2; i++ ) {

            if ( data[GROUP_DENS][i] == 0 ) continue;

            for ( j=0; j<GROUP_FILED_NBLOCKS; j++ ) {

                if ( !group_filed_present( j ) )
                    continue;
                if ( j==GROUP_DENS )
                    continue;
                data[j][i] /= data[GROUP_DENS][i];

            }

        }

        for ( i=0; i<PicSize2; i++ ) {

            if ( num[i] == 0 )
                continue;
            data[GROUP_DENS][i] /= num[i];
            data[GROUP_DENS][i] /= Time3;
            data[GROUP_DENS][i] /= RhoBaryon;
            //data[GROUP_DENS][i]  *= (( g2c.g ) / CUBE( g2c.cm ) / CUBE( Time ));

        }


        for ( i=0; i<6; i++ )
            img_props(i) = mass[i];
        img_props( 6 ) = r200;
        img_xmin =  -L;
        img_xmax =  L;
        img_ymin =  -L;
        img_ymax =  L;


        for( i=0; i<GROUP_FILED_NBLOCKS; i++ ) {

            if ( !group_filed_present( i ) )
                continue;

            if ( i == GROUP_RAD )
                continue;

            image.img = data[i];
            //conv( dL );
            get_group_filed_name( i, buf1 );
            make_group_output_filename( buf, buf1, g_index );
            write_img( buf, buf1 );

        }

        if ( group_filed_present( GROUP_RAD ) )
            if ( Redshift > 1e-10 ) {
                for ( i=0; i<SQR(PicSize); i++ ) {
                    data[GROUP_RAD][i]  = data[GROUP_RAD][i] / ( 4 * PI * SQR( LumDis*g2c.cm ) )
                                       //   / SQR( dL / ComDis  ) ;
                                       / SQR( dL / ComDis  / PI * 180 * 60 ) ;
                }

                img_xmin =  -L / ComDis / PI * 180 * 60;
                img_xmax =   L / ComDis / PI * 180 * 60;
                img_ymin =  -L / ComDis / PI * 180 * 60;
                img_ymax =   L / ComDis / PI * 180 * 60;
                image.img = data[GROUP_RAD];
                get_group_filed_name( GROUP_RAD, buf1 );
                make_group_output_filename( buf, buf1, g_index );
                write_img( buf, buf1 );

        }

    } // for index

    myfree( num );

    for( i=0; i<GROUP_FILED_NBLOCKS; i++ ) {

        if ( !group_filed_present( i ) )
            continue;
        myfree( data[i] );

    }

    put_sep0;

    reset_img();

}

void group_analysis() {

    int proj_tmp[3], i;
    char Sproj_tmp;

    for( i=0; i<3; i++ )
        proj_tmp[i] = Proj[i];
    Sproj_tmp = Sproj;

    for( i=0; i<3; i++ ) {
        proj_k = i;
        proj_i = ( i + 1  ) % 3;
        proj_j = ( i + 2  ) % 3;
        Sproj = i + 'x';
        group_plot();
    }

    for( i=0; i<3; i++ )
        Proj[i] = proj_tmp[i];
    Sproj = Sproj_tmp;

    if ( All.GroupEleSpec )
        group_electron_spectrum();

    put_sep0;

    if ( All.GroupSpec ) {
        group_spectrum();
        //group_spectrum_index();
    }

    if ( All.GroupTempProfile )
        group_temp_profile();

    if ( All.GroupTempStack )
        group_temp_stack();

    if ( All.GroupGasRatio )
        group_gas_ratio();

    reset_img();
    put_sep0;

}

