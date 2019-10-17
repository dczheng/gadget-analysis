#include "allvars.h"


void group_temp_profile() {
#ifdef GROUPTEMPPROFILE

    long p;
    struct group_properties g;
    int g_index, i, RN, *num, ii, ngbnum, x, y, z, k;
    double  *T, *data, Rmin, Rmax, dlogR, r;
    char buf[100], *dn="TempProfile";
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
    for ( g_index=0; g_index<Ngroup; g_index++ ) {

        if ( !group_present( g_index ) )
            break;

        g = Gprops[g_index];
        ngbnum = ngb( g.cm, Rmax*0.9, 0 );
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
            write_img( buf );
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

#endif
}

void group_temp_stack() {
#ifdef GROUPTEMPSTACK

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

    for ( g_index=0; g_index<Ngroup; g_index++ ) {
        g = Gprops[g_index];
        if ( g.mass < m[0] )
            continue;

        mi = 0;
        while( mi < MS-1 && g.mass > m[mi+1] ) mi ++;

        ngbnum = ngb( g.cm, All. GroupTempStackRmax * 1.1, 0 );
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

    for ( g_index=0; g_index<Ngroup; g_index++ ) {
        g = Gprops[g_index];
        if ( g.mass < m[0] )
            continue;
        mi = 0;
        while( mi < MS-1 && g.mass > m[mi+1] ) mi ++;

        ngbnum = ngb( g.cm, All. GroupTempStackRmax, 0 );
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

#endif
}

