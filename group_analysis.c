#include "allvars.h"

#define make_group_output_filename( buf, nstr, group_index ) \
    sprintf( buf, "%s/%s/%s_%.2f_%04i_%c.dat",\
            All.GroupDir, nstr, nstr, All.RedShift, group_index, All.Sproj );

int group_present( long index ) {

    if ( ( index >= All.GroupIndexMin ) &&
         ( index <= All.GroupIndexMax ) &&
         ( Gprops[index].mass >= All.GroupMassMin) )
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

inline double get_group_size( struct group_properties *g ) {
    return vmax( g->size[All.proj_i], g->size[All.proj_j] );
}

void group_flux( int nu_index, long index, double *flux, double *flux_nosr ) {

    double com_dis, lum_dis, L;
    struct group_properties *g;
    double size;

    g = &Gprops[index];
    com_dis = comoving_distance( All.Time );
    lum_dis = luminosity_distance( All.Time ) * g2c.cm;

    L = group_luminosity( nu_index, index );

    size = get_group_size( g );
    *flux_nosr = L / ( 4.0 * PI * SQR( lum_dis ) );
    *flux = *flux_nosr / ( SQR( size * 2  / com_dis ) );

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
    char buf[100];
    FILE *fd;

    qn = All.QNum;
    qmin = All.QMin;
    qmax = All.QMax;

    dlogq = log( qmax/qmin ) / ( qn-1 );

    sprintf( buf, "%s/ElectronSpectrum", All.GroupDir );
    create_dir( buf );

    sprintf( buf, "%s/ElectronSpectrum/EleSpec_%.2f.dat",
            All.GroupDir, All.RedShift );

    fd = fopen( buf, "w" );

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
        memset( f, 0, sizeof(double) * qn );
        while( p >= 0 ) {
            if ( P[p].Type == 0 ) {
                for( i=0; i<qn; i++ ) {
                    q = log( qmin ) + i * dlogq;
                    q = exp(q);
                    f[i] += particle_f( &SphP[p], q ) * SphP[p].Density;
                    rho += SphP[p].Density;
                }
                /*
                if ( SphP[p].CRE_qmax > 1e4 )
                    printf( "[%i] %g, %g, %g\n",
                            index,
                            SphP[p].CRE_Alpha,
                            SphP[p].CRE_qmin,
                            SphP[p].CRE_qmax);
                            */
                vmax2( qmax_max, SphP[p].CRE_qmax );
                vmin2( qmin_min, SphP[p].CRE_qmin, 1 );
            }
            p = FoFNext[p];
        }

        printf( "[%i] qmax_max: %g, qmin_min: %g\n", index, qmax_max, qmin_min );
        fprintf( fd, "%i  %g  ", index, Gprops[index].mass );
        for ( i=0; i<qn; i++ ) {
            fprintf( fd, "%g  ", f[i]/rho );
        }
        fprintf( fd, "\n" );
    }

    myfree( f );

    fclose( fd );

}

void group_spectrum() {

    int vN, i, index;
    double v, *flux, vmin, vmax, dv, *flux_nosr;
    char buf[100];
    FILE *fd1, *fd2;

    writelog( "group spectrum ...\n" );

    vN = All.NuNum;
    vmin = All.NuMin;
    vmax = All.NuMax;

    dv = log( vmax/vmin ) / ( vN-1 );

    sprintf( buf, "%s/Spectrum", All.GroupDir );
    create_dir( buf );

    sprintf( buf, "%s/Spectrum/Spec_%.2f.dat",
            All.GroupDir, All.RedShift );

    fd1 = fopen( buf, "w" );

    sprintf( buf, "%s/Spectrum/Spec_%.2f_nosr.dat",
            All.GroupDir, All.RedShift );

    mymalloc1( flux, sizeof(double) * vN );
    mymalloc1( flux_nosr, sizeof(double) * vN );

    fd2 = fopen( buf, "w" );

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

        for ( i=0; i<vN; i++ )
            group_flux( i, index, flux+i, flux_nosr+i );

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

    writelog( "group spectrum ... done.\n" );

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
    x = All.proj_i;
    y = All.proj_j;
    xo = PicS / 2;
    yo = PicS / 2;

    dv = log10( vmax/vmin) / vN;

    sprintf( buf, "%s/%s", All.GroupDir, spec_index_str );
    create_dir( buf );

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
        write_img1( buf, spec_index_str );

        image.img = spec_index_err;
        sprintf( buf, "%s/%s/%s_%.2f_%04i_%c.dat",
            All.GroupDir, spec_index_str, spec_index_err_str,
            All.RedShift, index, All.Sproj );
        write_img1( buf, spec_index_err_str );

    }

    myfree( v );
    myfree( spec );
    myfree( spec_index );
    myfree( spec_index_err );
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
        case GROUP_RADP:
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
        case GROUP_RADP:
            strcpy( buf, "RadioP" );
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

void output_group_info(  long index ) {

    long p;
    struct group_properties *g;
    FILE *fd;
    char buf[20];

    sprintf( buf, "group_%04li.dat", index );

    fd = fopen( buf, "w" );
    g = &Gprops[index];

    p = g->Head;

    while( p >= 0 ) {

        //printf( "%li\n", p );
        if ( P[p].Type == 0 ) {
            fprintf( fd, "%g %g %g  %g %g %g %g %g %g %g\n",
                //P[p].Pos[0],
                //P[p].Pos[1],
                //P[p].Pos[2],
                SphP[p].CR_C0 * pow( SphP[p].Density, (2.75-1)*0.3333 )
                * SphP[p].Density / guc.m_p
                / CUBE( g2c.cm ),

                SphP[p].CR_n0
                * SphP[p].Density / guc.m_p
                / CUBE( g2c.cm ),

                SphP[p].CR_E0
                * SphP[p].Density
                * g2c.erg / CUBE( g2c.cm ),

                get_B( p ),

                SphP[p].CRE_C * SphP[p].Density / guc.m_e
                / CUBE( g2c.cm ),

                SphP[p].CRE_Alpha,
                SphP[p].CRE_qmin,
                SphP[p].CRE_qmax,

                SphP[p].CRE_n * SphP[p].Density / guc.m_e
                / CUBE( g2c.cm ),

                SphP[p].CRE_e * SphP[p].Density
                * g2c.erg / CUBE( g2c.cm )
                );
        }

        p = FoFNext[p];

    }

    fclose( fd );

    //endrun( 20190123 );

}

void group_analysis() {

    long *npart, p;
    struct group_properties g;
    double L, dL, *mass, B, PP, nu, n, e,
            lum_dis, com_dis;// com_dis, lum_dis, h, ang;
    int *num, g_index, i, j, PicSize, x, y,
         xo, yo, PicSize2, pic_index, ii, jj;

    char buf[100], buf1[100];
    double *data[GROUP_FILED_NBLOCKS];

    if ( ThisTask_Local != 0 )
        return;

    if ( ThisTask == 0 )
        output_group_info( 0 );

    PicSize = All.PicSize;
    PicSize2 = All.PicSize2;
    x = All.proj_i;
    y = All.proj_j;
    xo = PicSize / 2;
    yo = PicSize / 2;

    mymalloc2( num, PicSize2 * sizeof( int ) );

    for( i=0; i<GROUP_FILED_NBLOCKS; i++ ) {

        if ( !group_filed_present( i ) )
            continue;

        mymalloc1( data[i], PicSize2 * sizeof( double ) );
        get_group_filed_name( i, buf1 );
        sprintf( buf, "%s/%s", All.GroupDir, buf1 );
        create_dir( buf );

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

        writelog( "analysis group: %i ...\n", g_index );

        g = Gprops[g_index];
        writelog( "center of mass: %g %g %g\n",
                g.cm[0], g.cm[1], g.cm[2] );

        p = g.Head;
        mass = g.mass_table;
        npart = g.npart;
        L = get_group_size( &g ) * 1.1;

        dL = 2 * L / PicSize;
        writelog( "npart: " );
        for ( i=0; i<6; i++ )
            writelog( "%li ", npart[i] );
        writelog( "\n" );

        writelog( "mass: " );
        for ( i=0; i<6; i++ ) {
            mass[i] *= 1e10;
            writelog( "%g ", mass[i] );
        }

        writelog( "\n" );
        writelog( "L: %g, dL:%g\n", 2*L, dL );

        p = g.Head;

        for ( i=0; i<g.Len; i++, p=FoFNext[p] ) {
            if ( P[p].Type != 0 )
                continue;

            /*
            printf( "[%li] Pos: ( %g, %g, %g )\n",
                    p, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2] );
            continue;
            */

            ii = PERIODIC( P[p].Pos[x] - g.cm[x] ) / dL + xo;
            jj = PERIODIC( P[p].Pos[y] - g.cm[y] ) / dL + yo;

            /*
            if ( ii < 0 || ii >= PicSize ||
                    jj < 0 || jj >= PicSize )
                    printf( "%i, %i\n", ii, jj );
                    */

            //continue;

            check_picture_index( ii );
            check_picture_index( jj );
            pic_index= ii*PicSize + jj;

            num[ pic_index ] ++;

            if ( group_filed_present( GROUP_DENS ) )
                data[GROUP_DENS][pic_index] += SphP[p].Density;

            if ( group_filed_present( GROUP_TEMP ) )
                data[GROUP_TEMP][pic_index] += SphP[p].Temp * SphP[p].Density;

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
                e = SphP[p].CRE_e * SphP[p].Density;
                e /= g2c.erg / CUBE( g2c.cm );
                data[GROUP_CREE][pic_index] += e * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_CREQMIN ) ) {
                data[GROUP_CREQMIN][pic_index] += SphP[p].CRE_qmin * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_CREQMAX ) ) {
                data[GROUP_CREQMAX][pic_index] += SphP[p].CRE_qmax * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_RAD ) ) {
                nu = All.Freq * 1e6;
                PP = particle_radio( nu, p );
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
                if ( j==GROUP_RADP )
                    continue;

                data[j][i] /= data[GROUP_DENS][i];

            }

        }

        for ( i=0; i<PicSize2; i++ ) {

            if ( num[i] == 0 )
                continue;
            data[GROUP_DENS][i] /= num[i];
            data[GROUP_DENS][i] /= CUBE( All.Time );
            data[GROUP_DENS][i] /= All.RhoBaryon;
            //data[GROUP_DENS][i]  *= (( g2c.g ) / CUBE( g2c.cm ) / CUBE( All.Time ));

        }


        for ( i=0; i<6; i++ )
            img_props(i) = mass[i];
        img_xmin =  -L;
        img_xmax =  L;
        img_ymin =  -L;
        img_ymax =  L;


        for( i=0; i<GROUP_FILED_NBLOCKS; i++ ) {

            if ( !group_filed_present( i ) )
                continue;
            if( i== GROUP_RADP) continue;

            image.img = data[i];
            //conv( dL );
            get_group_filed_name( i, buf1 );
            make_group_output_filename( buf, buf1, g_index );
            write_img1( buf, buf1 );

        }

        if ( group_filed_present( GROUP_RADP ) )
            if ( All.RedShift > 1e-10 ) {
                lum_dis = luminosity_distance( All.Time ) * ( g2c.cm );
                com_dis = comoving_distance( All.Time );
     //           printf( "z:%g lum_dis: %g\n", All.RedShift, lum_dis);
                /*
                ang_dis = angular_distance( All.Time );
                h = All.SofteningTable[0];
                ang = h / com_dis / PI * 180;
                */
                for ( i=0; i<SQR(PicSize); i++ ) {
                    data[GROUP_RADP][i]  = data[GROUP_RAD][i] / ( 4 * PI * SQR( lum_dis ) ) * 1e23 * CUBE( KPC ) * 1e3;
                }

                img_xmin =  -L / com_dis / PI * 180 * 60;
                img_xmax =  L / com_dis / PI * 180 * 60;
                img_ymin =  -L / com_dis / PI * 180 * 60;
                img_ymax =  L / com_dis / PI * 180 * 60;
                image.img = data[GROUP_RADP];
                get_group_filed_name( GROUP_RADP, buf1 );
                make_group_output_filename( buf, buf1, g_index );
                write_img1( buf, buf1 );

        }

    } // for index

    myfree( num );

    for( i=0; i<GROUP_FILED_NBLOCKS; i++ ) {

        if ( !group_filed_present( i ) )
            continue;
        myfree( data[i] );

    }

    put_sep0;

    if ( All.GroupEleSpec )
        group_electron_spectrum();

    put_sep0;

    if ( All.GroupSpec ) {
        group_spectrum();
        //group_spectrum_index();
    }

    put_sep0;

}

