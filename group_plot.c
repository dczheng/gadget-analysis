#include "allvars.h"

int group_filed_present( enum group_fields blk ) {

    switch( blk ) {
        case GROUP_DENS:
                return 1;
        case GROUP_TEMP:
            if ( All.GroupTemp )
                return 1;
            return 0;
        case GROUP_U:
            if ( All.GroupU )
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
        case GROUP_U:
            strcpy( buf, "InternalEnergy" );
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

void group_plot() {

    long p;
    struct group_properties g;
    double L, dL, *mass, PP, n, r200, dlognu, nu_x, w; 
    int *num, g_index, i, j, x, y,
         xo, yo, pic_index, ii, jj, nu_i;

    char buf[100], buf1[100];
    double *data[GROUP_FILED_NBLOCKS];
//#define GROUP_PLOT_DEBUG

#ifdef GROUP_PLOT_DEBUG
    FILE *fd_d;
    int ti, tj;
    fd_d = myfopen( "w", "group_plot-test.dat" );
#endif

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
        while( p>=0 ){
            if ( P[p].Type != 0 ) {
                p = FoFNext[p];
                continue;
            }

            //printf( "[%li] Pos: ( %g, %g, %g )\n",
            //        p, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2] );
            //continue;

            ii = PERIODIC_HALF( P[p].Pos[x] - g.cm[x] ) / dL + xo;
            jj = PERIODIC_HALF( P[p].Pos[y] - g.cm[y] ) / dL + yo;

/*
            if ( get_B(p) * 1e6 > 100 ) {
                p = FoFNext[p];
                continue;
            }
*/
#ifdef GROUP_PLOT_DEBUG
            //for( ti=-1; ti<2; ti++ )
            //    for( tj=-1; tj<2; tj++ ) {
             //       if ( ii==xo+ti && jj==yo+tj )
                        fprintf( fd_d, "%g\n", get_B(p)*1e6 );
            //}
            p = FoFNext[p];
            continue;
#endif

            //if ( ii < 0 || ii >= PicSize ||
            //        jj < 0 || jj >= PicSize )
            //        printf( "%i, %i\n", ii, jj );
            //
            //continue;

            check_picture_index( ii );
            check_picture_index( jj );
            pic_index= ii*PicSize + jj;

            num[pic_index] ++;

            if ( All.GroupDensityWeighted )
                w = SphP[p].Density;
            else
                w = 1;


            if ( group_filed_present( GROUP_DENS ) )
                data[GROUP_DENS][pic_index] += SphP[p].Density;

            if ( group_filed_present( GROUP_TEMP ) )
                data[GROUP_TEMP][pic_index] += SphP[p].Temp * w;

            if ( group_filed_present( GROUP_SFR ) )
                data[GROUP_SFR][pic_index] += SphP[p].sfr * w;

            if ( group_filed_present( GROUP_U ) )
                data[GROUP_U][pic_index] += SphP[p].u * w;

            if ( group_filed_present( GROUP_MAG ) )
                data[GROUP_MAG][pic_index] += get_B(p) * 1e6 * w;

            if ( group_filed_present( GROUP_MACH ) )
                data[GROUP_MACH][pic_index] += SphP[p].MachNumber * w;

            if ( group_filed_present( GROUP_CREN ) ) {
                n = SphP[p].CRE_n *  SphP[p].Density / guc.m_e;
                n /= CUBE( g2c.cm );
                data[GROUP_CREN][pic_index] += n * w;
            }

            if ( group_filed_present( GROUP_CREALPHA ) )
                data[GROUP_CREALPHA][pic_index] += SphP[p].CRE_Alpha * w;

            if ( group_filed_present( GROUP_CREC ) )
                data[GROUP_CREC][pic_index] +=  SphP[p].CRE_C * w;

            if ( group_filed_present( GROUP_CREE ) ) {
                /*
                e = SphP[p].CRE_e * SphP[p].Density;
                e /= g2c.erg / CUBE( g2c.cm );
                data[GROUP_CREE][pic_index] += e * SphP[p].Density;
                */
                data[GROUP_CREE][pic_index] += SphP[p].CRE_e / SphP[p].u * w;
            }

            if ( group_filed_present( GROUP_CREQMIN ) )
                data[GROUP_CREQMIN][pic_index] += SphP[p].CRE_qmin * w;

            if ( group_filed_present( GROUP_CREQMAX ) )
                data[GROUP_CREQMAX][pic_index] += SphP[p].CRE_qmax * w;

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
                data[GROUP_RAD][pic_index] += PP * w;
            }

            p = FoFNext[p];
        } // while
#ifdef GROUP_PLOT_DEBUG
        break;
#endif

        for ( i=0; i<PicSize2; i++ ) {

            if ( data[GROUP_DENS][i] == 0 ) continue;

            for ( j=0; j<GROUP_FILED_NBLOCKS; j++ ) {

                if ( !group_filed_present( j ) )
                    continue;
                if ( j==GROUP_DENS )
                    continue;

                if ( All.GroupDensityWeighted )
                    data[j][i] /= data[GROUP_DENS][i];
                else
                    data[j][i] /= num[i];

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
                for ( i=0; i<PicSize2; i++ ) {
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
#ifdef GROUP_PLOT_DEBUG
        //break;
#endif

    } // for index

    myfree( num );

    for( i=0; i<GROUP_FILED_NBLOCKS; i++ ) {

        if ( !group_filed_present( i ) )
            continue;
        myfree( data[i] );

    }

#ifdef GROUP_PLOT_DEBUG
        fclose( fd_d );
        endruns( "group_plot-test" );
#endif

    reset_img();

}

