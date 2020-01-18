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
            return SphP[p].Density / Time3 / RhoBaryon;
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
                return log(f/f1) / log(All.GroupRadFreq1/All.GroupRadFreq);
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

#define ALTINDEX
#if defined(GROUPTEMP) || defined(GROUPU) || defined(GROUPSFR) || defined(GROUPB) || defined(GROUPMACH) || defined(GROUPCRE) || defined(GROUPRAD) || defined(GROUPRADINDEX) || defined(GROUPDENSITY)
    long p, *pidx, nidx;
    struct group_properties g;
    double L, dL, *mass, r200, area, fac_to_arcmin, Lz, px, py, h, fac,
            u, r, m, rho, wk, hinv, hinv3, hinv4, dwk; 
    int  g_index, i, j, x, y, k,
         pic_index, z, proj_tmp, i1, i2, j1, j2;

    char buf[100], buf1[100];
    double *data[GROUP_FIELD_NBLOCKS];
    put_header( "group_plot" );

    fac_to_arcmin = 1 / PI * 180 * 60;

   mymalloc1( pidx, sizeof(long)*N_Gas );
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

    for( i=0; i<GROUP_FIELD_NBLOCKS; i++ ) {

        if ( !group_filed_present( i ) )
            continue;

        mymalloc1( data[i], PicSize2 * sizeof( double ) );
        get_group_filed_name( i, buf1 );
        create_dir( "%s%s", GroupDir, buf1 );

    }
#ifdef GROUPDENSITYWEIGHTED
    double *ww;
    mymalloc1( ww, PicSize2 * sizeof(double) );
#endif
    
    proj_tmp = All.ProjectDirection;
    for( z=0; z<3; z++ ) {

        x = (z+1) % 3;
        y = (z+2) % 3;

        writelog( "%c (%c %c)\n", 'x'+z, 'x'+x, 'x'+y );
       
        Sproj = 'x'+z;
        All.ProjectDirection = z;
        for ( g_index=0; g_index<Ngroup; g_index++ ) {
    
            if ( !group_present( g_index ) )
                break;

            if ( (g_index * 3 + z) % NTask_Local != ThisTask_Local )
                continue;
    
            for( i=0; i<GROUP_FIELD_NBLOCKS; i++ ) {
    
                if ( !group_filed_present(i) )
                    continue;
                memset( data[i], 0, PicSize2 * sizeof( double ) );
            }
#ifdef GROUPDENSITYWEIGHTED
            memset( ww, 0, PicSize2 * sizeof( double ) );
#endif

            g = Gprops[g_index];
            p = g.Head;
            mass = g.mass_table;
            r200 = g.vr200;
    
#ifndef GROUPFIXEDSIZE
            L = vmax( g.size[x], g.size[y] );
            L *= 1.1;
#else
            L = vmax(L3[x], L3[y]);
#endif
    
            dL = 2 * L / PicSize;
            p = g.Head;
            Lz = All.GroupSizeZ;
            if ( Lz == 0 )
                Lz = g.size[z];

            for( p=0, nidx=0; p<N_Gas; p++ ) {
                if (
                    ( NGB_PERIODIC(P[p].Pos[x] - g.cm[x]) > L ) ||
                    ( NGB_PERIODIC(P[p].Pos[y] - g.cm[y]) > L ) ||
                      NGB_PERIODIC(P[p].Pos[z] - g.cm[z]) > Lz ) 
                        continue;
                pidx[nidx++] = p;
            }
            
            for ( p=0; p<nidx; p++ ) {
            
                h = SphP[pidx[p]].Hsml;                             
                m = P[pidx[p]].Mass;       
                rho = SphP[pidx[p]].Density;    
                px = PERIODIC_HALF(P[pidx[p]].Pos[x] - g.cm[x]) + L;
                py = PERIODIC_HALF(P[pidx[p]].Pos[y] - g.cm[y]) + L;
            
                j1 = (px-h)/dL-1;
                j2 = (px+h)/dL+1;
                i1 = (py-h)/dL-1;
                i2 = (py+h)/dL+1;
                //printf( "%i, %i, %i, %i\n", i1, i2, j1, j2 );
            
                for ( i=i1; i<=i2; i++ ) {
                    for (j=j1; j<=j2; j++) {
                        pic_index = j * PicSize + i;
                        if ( pic_index < 0 || pic_index > PicSize*PicSize-1 )
                            continue;
                        r = sqrt( pow( i*dL - dL/2 - py, 2 ) + 
                                  pow( j*dL - dL/2 - px, 2 ) );
                        //printf( "%g\n", r );
                        if ( r > h )
                            continue;
            
                        u = r / h;
                        kernel_hinv( h, &hinv, &hinv3, &hinv4 );
                        kernel_main( u, hinv3, hinv4, &wk, &dwk, -1 );
                        fac = wk * m / rho;

                        for( k=0; k<GROUP_FIELD_NBLOCKS; k++ ) {
                            if ( !group_filed_present(k) )
                                continue;
#ifdef ALTINDEX
                            if ( k == GROUP_RADINDEX )
                                continue;
#endif
                            data[k][pic_index] += get_group_filed_data(k, pidx[p]) * fac;
                        }
                    }
                }
                
            }

#ifdef ALTINDEX
#ifdef GROUPRADINDEX
            for( i=0; i<PicSize2; i++ ) {
                if ( data[GROUP_RAD1][i] > 0 ) {
                    data[GROUP_RADINDEX][i] =
                    log( data[GROUP_RAD][i] / data[GROUP_RAD1][i]) /
                                            log( All.GroupRadFreq1/All.GroupRadFreq );
                }
            }
#endif
#endif
            for( i=0; i<GROUP_FIELD_NBLOCKS; i++ ) {
    
                if ( !group_filed_present( i ) )
                    continue;
    
                if ( i == GROUP_RAD || i == GROUP_RAD1 ) {
    
                    area = get_solid_angle( dL, dL, ComDis ) * SQR(fac_to_arcmin);
                    if ( Redshift > 1e-10 ) {
                        for ( j=0; j<PicSize2; j++ ) {
                            data[i][j]  = data[i][j] / ( 4 * PI * SQR( LumDis*g2c.cm ) )
                                           / area;
                        }
    
                        img_xmin =  img_ymin = -L / ComDis * fac_to_arcmin;
                        img_xmax =  img_ymax = -img_xmin;
                        img_props(7) = r200 / ComDis * fac_to_arcmin;
                        image.img = data[i];
                        get_group_filed_name( i, buf1 );
                        make_group_output_filename( buf, buf1, g_index );
                        write_img( buf );
    
                        continue;
                    }
                }

                for ( j=0; j<6; j++ )
                    img_props(j) = mass[j];
                img_props(7) = r200;
                img_xmin =  -L;
                img_xmax =  L;
                img_ymin =  -L;
                img_ymax =  L;
    
                image.img = data[i];
                //conv( dL );
                get_group_filed_name( i, buf1 );
                make_group_output_filename( buf, buf1, g_index );
                write_img( buf );
            }
        break;
        } // for gindex

    } // for z

    All.ProjectDirection = proj_tmp;
    Sproj = proj_tmp + 'x';
    
#ifdef GROUPDENSITYWEIGHTED
    myfree( ww );
#endif

    for( i=0; i<GROUP_FIELD_NBLOCKS; i++ ) {
        if ( !group_filed_present( i ) )
            continue;
        myfree( data[i] );
    }

    myfree( pidx );

    reset_img();
    put_end();

#endif
}


