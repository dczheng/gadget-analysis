#include "allvars.h"


int group_filed_present( enum group_fields blk ) {

    switch( blk ) {
        case GROUP_DENS:
                return 1;
        case GROUP_TEMP:
#ifdef GROUPTEMP
                return 1;
#endif
            return 0;
        case GROUP_U:
#ifdef GROUPU
                return 1;
#endif
            return 0;
        case GROUP_SFR:
#ifdef GROUPSFR
                return 1;
#endif
            return 0;
        case GROUP_MAG:
#ifdef GROUPB
                return 1;
#endif
            return 0;
        case GROUP_MACH:
#ifdef GROUPMACH
                return 1;
#endif
            return 0;
        case GROUP_CREN:
        case GROUP_CREC:
        case GROUP_CREE:
        case GROUP_CREALPHA:
        case GROUP_CREQMIN:
        case GROUP_CREQMAX:
#ifdef GROUPCRE
                return 1;
#endif
            return 0;
        case GROUP_RAD:
        case GROUP_RAD1:
#ifdef GROUPRAD
                return 1;
#endif
        case GROUP_RADINDEX:
#ifdef GROUPRADINDEX
                return 1;
#endif
            return 0;
        default:
            return 0;
    }
}

double get_group_filed_data( enum group_fields blk, long p ) {

    switch( blk ) {
        case GROUP_DENS:
            return SphP[p].Density;
#ifdef GROUPTEMP
        case GROUP_TEMP:
            return SphP[p].Temp;
#endif
#ifdef GROUPU
        case GROUP_U:
            return SphP[p].u;
#endif
#ifdef GROUPSFR
        case GROUP_SFR:
            return SphP[p].sfr;
#endif
#ifdef GROUPB
        case GROUP_MAG:
            return get_B(p) * 1e6;
#endif
#ifdef GROUPMACH
        case GROUP_MACH:
            return SphP[p].MachNumber;
#endif
#ifdef GROUPCRE
        case GROUP_CREN:
            return SphP[p].CRE_n *  SphP[p].Density / guc.m_e / CUBE( g2c.cm );
        case GROUP_CREC:
            return SphP[p].CRE_C;
        case GROUP_CREE:
            return SphP[p].CRE_e / SphP[p].u;
        case GROUP_CREALPHA:
            return SphP[p].CRE_Alpha;
        case GROUP_CREQMIN:
            return SphP[p].CRE_qmin;
        case GROUP_CREQMAX:
            return SphP[p].CRE_qmax;
#endif
#ifdef GROUPRAD
        case GROUP_RAD:
             return get_particle_radio( p, All.GroupRadFreq );
        case GROUP_RAD1:
             return get_particle_radio( p, All.GroupRadFreq1 );
#endif
#ifdef GROUPRAD
        double f, f1;
        case GROUP_RADINDEX:
            f = get_particle_radio( p, All.GroupRadFreq );
            f1 = get_particle_radio( p, All.GroupRadFreq1 );
            if ( f1 != 0 )
                return log(f/f1) /
                    log(All.GroupRadFreq1/All.GroupRadFreq);
            return 0;
#endif
       default:
            endruns( "can't occur !!!" );
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
        case GROUP_RAD1:
            strcpy( buf, "Radio1" );
            break;
        case GROUP_RADINDEX:
            strcpy( buf, "RadioIndex" );
            break;

        default:
            endruns( "can't occur !!!" );
    }
}

void group_plot() {

#if defined(GROUPTEMP) || defined(GROUPU) || defined(GROUPSFR) || defined(GROUPB) || defined(GROUPMACH) || defined(GROUPCRE) || defined(GROUPRAD) || defined(GROUPRADINDEX)
    long p;
    struct group_properties g;
    double L, dL, *mass, r200, w, fac_arc; 
    int *num, g_index, i, j, x, y,
         xo, yo, pic_index, ii, jj, z;

    char buf[100], buf1[100];
    double *data[GROUP_FIELD_NBLOCKS];
    put_header( "group_plot" );

//#define GROUP_PLOT_DEBUG
#ifdef GROUP_PLOT_DEBUG
    FILE *fd_d;
    int ti, tj;
    fd_d = myfopen( "w", "group_plot-test.dat" );
#endif

    fac_arc = 180 * 60 * 60 / PI;

#ifdef GROUPFIXEDSIZE
    double L3[3];
    if ( All.GroupSize > 0 ) {
        for( i=0; i<3; i++ )
            L3[i] = All.GroupSize;
    }
    else{
    
        for( z=0; z<3; z++ ) L3[z] = 0;
        for ( g_index=0; g_index<Ngroup; g_index++ ) {
            if ( !group_present( g_index ) )
                break;
            g = Gprops[g_index];
    
            for( z=0; z<3; z++ ) {
                vmax2( L3[z], g.size[(z+1)%3]);
                vmax2( L3[z], g.size[(z+2)%3]);
            }
        }
        for( z=0; z<3; z++ ) L3[z] *= 1.1;
    }
    writelog( "L: %g %g %g\n", L3[0], L3[1], L3[2] );
#endif

    mymalloc2( num, PicSize2 * sizeof( int ) );

    for( i=0; i<GROUP_FIELD_NBLOCKS; i++ ) {

        if ( !group_filed_present( i ) )
            continue;

        mymalloc1( data[i], PicSize2 * sizeof( double ) );
        get_group_filed_name( i, buf1 );
        create_dir( "%s%s", GroupDir, buf1 );

    }

    xo = PicSize / 2;
    yo = PicSize / 2;
    
    for( z=0; z<3; z++ ) {

        x = (z+1) % 3;
        y = (z+2) % 3;

        writelog( "%c (%c %c)\n", 'x'+z, 'x'+x, 'x'+y );
       
        for ( g_index=0; g_index<Ngroup; g_index++ ) {
    
            if ( !group_present( g_index ) )
                break;

            if ( (g_index * 3 + z) % NTask_Local != ThisTask_Local )
                continue;
    
            memset( num, 0, PicSize2 * sizeof( int ) );
            for( i=0; i<GROUP_FIELD_NBLOCKS; i++ ) {
    
                if ( !group_filed_present(i) )
                    continue;
                memset( data[i], 0, PicSize2 * sizeof( double ) );
    
            }
//#define PRINT_GROUP_INFO
            g = Gprops[g_index];
            p = g.Head;
            mass = g.mass_table;
            r200 = g.vr200;
    
#ifndef GROUPFIXEDSIZE
           L = vmax( g.size[x], g.size[y] );
           L *= 1.1;
#else
           L = L3[z];
#endif
    
           dL = 2 * L / PicSize;
#ifdef PRINT_GROUP_INFO
            writelog( "group plot[%c]: %i ...\n", z+'x', g_index );
            writelog( "center of mass: %g %g %g\n",
                    g.cm[0], g.cm[1], g.cm[2] );
    
            writelog( "r200: %g\n", r200 );
            writelog( "npart: " );
            for ( i=0; i<6; i++ )
                writelog( "%li ", g.npart[i] );
            writelog( "\n" );
    
            writelog( "mass: " );
            for ( i=0; i<6; i++ ) {
                writelog( "%g ", mass[i] );
            }
            writelog( "\n" );
            writelog( "L: %g, dL:%g\n", 2*L, dL );
#endif
    
            p = g.Head;
            while( p>=0 ){
                if ( P[p].Type != 0 ) {
                    p = FoFNext[p];
                    continue;
                }
                ii = PERIODIC_HALF( P[p].Pos[x] - g.cm[x] ) / dL + xo;
                jj = PERIODIC_HALF( P[p].Pos[y] - g.cm[y] ) / dL + yo;
#ifdef GROUP_PLOT_DEBUG
                //for( ti=-1; ti<2; ti++ )
                //    for( tj=-1; tj<2; tj++ ) {
                 //       if ( ii==xo+ti && jj==yo+tj )
                            fprintf( fd_d, "%g %g %g\n", get_B(p)*1e6,
                            SphP[p].MachNumber, SphP[p].Hsml );
                //}
                //p = FoFNext[p];
                //continue;
#endif
                check_picture_index( ii );
                check_picture_index( jj );
                pic_index= ii*PicSize + jj;
    
                num[pic_index] ++;
    
#ifdef GROUPDENSITYWEIGHTED
                w = SphP[p].Density;
#else
                w = 1;
#endif
                for( i=0; i<GROUP_FIELD_NBLOCKS; i++ ) {
                    if ( !group_filed_present(i) )
                        continue;
                    if ( i == GROUP_DENS )
                        data[i][pic_index] += get_group_filed_data(i, p);
                    else
                        data[i][pic_index] += get_group_filed_data(i, p) * w;
#ifdef GROUP_PLOT_DEBUG
                get_group_filed_name( i, buf1 );
                printf( "%s: %g ", buf1, get_group_filed_data( i, p ) );
#endif
                }
#ifdef GROUP_PLOT_DEBUG
                printf( "\n" );
#endif
                p = FoFNext[p];
            } 
#ifdef GROUP_PLOT_DEBUG
            endruns( "group-plot-debug" );
#endif
    
            for ( i=0; i<PicSize2; i++ ) {
    
                if ( data[GROUP_DENS][i] == 0 ) continue;
    
                for ( j=0; j<GROUP_FIELD_NBLOCKS; j++ ) {
    
                    if ( !group_filed_present( j ) )
                        continue;
                    if ( j==GROUP_DENS )
                        continue;
    
#ifdef GROUPDENSITYWEIGHTED
                        data[j][i] /= data[GROUP_DENS][i];
#else
                        data[j][i] /= num[i];
#endif
    
                }
    
            }
    
            for ( i=0; i<PicSize2; i++ ) {
    
                if ( num[i] == 0 )
                    continue;
                data[GROUP_DENS][i] /= num[i];
                data[GROUP_DENS][i] /= Time3;
                data[GROUP_DENS][i] /= RhoBaryon;
                data[GROUP_DENS][i]  *= (( g2c.g ) / CUBE( g2c.cm ) / CUBE( Time ));
    
            }
    
    
            for ( i=0; i<6; i++ )
                img_props(i) = mass[i];
            img_props(7) = r200;
            img_xmin =  -L;
            img_xmax =  L;
            img_ymin =  -L;
            img_ymax =  L;
    
            for( i=0; i<GROUP_FIELD_NBLOCKS; i++ ) {
    
                if ( !group_filed_present( i ) )
                    continue;
    
                if ( i == GROUP_RAD || i == GROUP_RAD1 ) {
    
                    if ( Redshift > 1e-10 ) {
                        for ( j=0; j<PicSize2; j++ ) {
                            data[i][j]  = data[i][j] / ( 4 * PI * SQR( LumDis*g2c.cm ) )
                                           //   / SQR( dL / ComDis  ) ;
                                           / SQR( dL / ComDis  / PI * 180 * 60 ) ;
                        }
    
                        img_xmin =  img_ymin = -L / ComDis * fac_arc;
                        img_xmax =  img_ymax = -img_xmin;
                        img_props(7) = r200 / ComDis * fac_arc;
                        image.img = data[i];
                        get_group_filed_name( i, buf1 );
                        make_group_output_filename( buf, buf1, g_index, z+'x' );
                        write_img( buf );
    
                        continue;
                    }
                }
    
                image.img = data[i];
                //conv( dL );
                get_group_filed_name( i, buf1 );
                make_group_output_filename( buf, buf1, g_index, z+'x' );
                write_img( buf );
    
            }
    
#ifdef GROUP_PLOT_DEBUG
            //break;
#endif
    
        } // for index

    } // for z
    
    myfree( num );

    for( i=0; i<GROUP_FIELD_NBLOCKS; i++ ) {
        if ( !group_filed_present( i ) )
            continue;
        myfree( data[i] );
    }

#ifdef GROUP_PLOT_DEBUG
        fclose( fd_d );
        endruns( "group_plot-test" );
#endif

    reset_img();
    put_end();

#endif
}
