#include "allvars.h"

void slice_init() {

    long i;
    int flag;

    put_header( "determine slice info" );

    SliceStart[0] = All.SliceStartX;
    SliceStart[1] = All.SliceStartY;
    SliceStart[2] = All.SliceStartZ;
    SliceEnd[0] = All.SliceEndX;
    SliceEnd[1] = All.SliceEndY;
    SliceEnd[2] = All.SliceEndZ;

    for( i=0; i<3; i++ )
        if ( SliceEnd[i] == 0 ) {
            SliceStart[i] = 0;
            SliceEnd[i] = BoxSize;
        }

    writelog( "Slice: " );
    for( i=0; i<3; i++ )
        writelog( "(%g, %g) ", SliceStart[i], SliceEnd[i] );
    writelog( "\n" );

    flag = 0;
    for( i=0; i<3; i++ )
        if ( SliceStart[i] != 0 || SliceEnd[i] != BoxSize ) {
            flag = 1;
            break;
        }

   if ( !flag ) {
       NumPartInSlice = NumPart;
       N_GasInSlice = N_Gas;
   }
   else {
   NumPartInSlice = N_GasInSlice = 0;
   for( i=0; i<NumPart; i++ ) 
       if ( InSlice(i) ) {
           NumPartInSlice++;
           if ( P[i].Type == 0 )
               N_GasInSlice ++;
       }
   }

   if ( SliceEnd[proj_i] - SliceStart[proj_i] !=
        SliceEnd[proj_j] - SliceStart[proj_j])
       endruns( "invalid slice" );

   SliceL = SliceEnd[proj_i] - SliceStart[proj_i];

   writelog( "NumPartInSlice: %li\n", NumPartInSlice );
   writelog( "N_GasInSlice: %li\n", N_GasInSlice );
   writelog( "SliceL: %g\n", SliceL );
   put_end();

}


int slice_field_present( enum  group_fields blk ) {

    switch( blk ) {
        case GROUP_MAG:
#ifdef BSLICE
        return 1;
#endif
        return 0;

        case GROUP_DENS:
#ifdef DENSITYSLICE
        return 1;
#endif
        return 0;

        case GROUP_MACH:
#ifdef MACHSLICE
        return 1;
#endif
        return 0;

        case GROUP_TEMP:
#ifdef TEMPSLICE
        return 1;
#endif
        return 0;

        case GROUP_CREE:
#ifdef CREESLICE
        return 1;
#endif
        return 0;

        case GROUP_CREN:
#ifdef CRENSLICE
        return 1;
#endif
        return 0;
        case GROUP_RAD:
#ifdef RADSLICE
        return 1;
#endif
        return 0;

        default:
            return 0;
    }

}

double get_slice_field_data( enum group_fields blk, long p ) {

    switch( blk ) {
        case GROUP_DENS:
#ifdef DENSITYSLICE
        return SphP[p].Density / Time3 / RhoBaryon;
#endif
        case GROUP_MAG:
#ifdef BSLICE
        return get_B( p ) * 1e6;
#endif

        case GROUP_MACH:
#ifdef MACHSLICE
        return SphP[p].MachNumber;
#endif
        case GROUP_TEMP:
#ifdef TEMPSLICE
        return SphP[p].Temp;
#endif
        case GROUP_CREE:
#ifdef CREESLICE
        return SphP[p].CRE_e / SphP[p].u;
#endif
        case GROUP_CREN:
#ifdef CRENSLICE
        return SphP[p].CRE_n * SphP[p].Density / guc.m_e / CUBE( g2c.cm );
#endif
        case GROUP_RAD:
#ifdef RADSLICE
        return get_particle_radio(i, All.RadSliceFreq) * 
                1.0 / (4.0 * PI * SQR( LumDis * g2c.cm )) / ( SQR(SliceL/PicSize) / SQR(ComDis) );
#endif
        default:
            endruns( "can't occur !!!"  );
            return 0;
    }

}

void get_slice_field_name( enum group_fields blk, char *buf ) {

    switch( blk ) {
        case GROUP_DENS:
#ifdef DENSITYSLICE
        sprintf( buf, "Density" );
        break;
#endif
        case GROUP_MAG:
#ifdef BSLICE
        sprintf( buf, "MagneticField" );
        break;
#endif

        case GROUP_MACH:
#ifdef MACHSLICE
        sprintf( buf, "MachNumber" );
        break;
#endif
        case GROUP_TEMP:
#ifdef TEMPSLICE
        sprintf( buf, "Temperature" );
        break;
#endif
        case GROUP_CREE:
#ifdef CREESLICE
        sprintf( buf, "cre_e" );
        break;
#endif
        case GROUP_CREN:
#ifdef CRENSLICE
        sprintf( buf, "cre_n" );
        break;
#endif
        case GROUP_RAD:
#ifdef RADSLICE
        sprintf( buf, "radio_%.2f", All.RadSliceFreq );
        break;
#endif
        default:
            endruns( "can't occur !!!"  );
    }

}

void field_slice( double *data, char *name, long N, double *weight) {

    char buf[100];

    //sprintf( buf, "`%s` slice\n", name );
    //writelog( buf );

    create_dir( "%s%s", OutputDir, name );
    sprintf( buf, "%s%s/%s_%03i.dat", OutputDir, name, name, SnapIndex );

    reset_img();
    data_to_grid2d( data, image.img, N, PicSize, SliceL, weight );
    img_xmin = 0;
    img_xmax = SliceL;
    img_ymin = 0;
    img_ymax = SliceL;

    write_img( buf );

}

void slice() {
    long i, idx;
    int k, flag;
    double *data, *weight;
    char buf[100];

    flag = 0;

    for(k=0; k<GROUP_FIELD_NBLOCKS; k++) {
        if ( !slice_field_present( k ) )
            continue;
        flag = 1;
        break;
    }

    if ( !flag )
        return;

    mymalloc2( data, 3 * sizeof( double ) * N_GasInSlice );
    mymalloc2( weight, sizeof( double ) * N_GasInSlice );
    for( k=0; k<GROUP_FIELD_NBLOCKS; k++ ) {

        memset( data, 0, sizeof(double) * 3 * N_GasInSlice );
        if ( !slice_field_present( k ) )
            continue;

        idx = 0;
        for ( i=0; i<N_Gas; i++ ) {
            if ( !InSlice(i) )
                continue;
            data[3*idx] = P[i].Pos[proj_i] - SliceStart[proj_i];
            data[3*idx+1] = P[i].Pos[proj_j] - SliceStart[proj_j];
            data[3*idx+2] = get_slice_field_data( k, i );
            weight[idx] = SphP[i].Density;
            idx++;
        }
        get_slice_field_name( k, buf );
        /*
        for( i=0; i<10; i++ ) {
            printf( "%g %g %g\n", data[i*3], data[i*3+1], data[i*3+2] );
        }
        endruns("");
        */
        field_slice( data, buf, N_GasInSlice, NULL );
        /*
        if ( k==GROUP_DENS )
            field_slice( data, buf, N_GasInSlice, NULL );
        else
            field_slice( data, buf, N_GasInSlice, weight );
        */
    }

    myfree( data );
    myfree( weight );
}

void density_slice() {
#ifdef DENSITYSLICE
    int i;
    long index;
    double *data;

    mymalloc2( data, sizeof( double ) * N_GasInSlice * 3 );

    index = 0;
    for ( i=0; i<N_Gas; i++ ) {
        if ( !InSlice(i) )
            continue;
        if ( SphP[i].Temp >= 1e7 ) {
            data[3*index] = P[i].Pos[proj_i] - SliceStart[proj_i];
            data[3*index+1] = P[i].Pos[proj_j] - SliceStart[proj_j];
            data[3*index+2] = SphP[i].Density / Time3 / RhoBaryon;
            index ++;
        }
    }
    field_slice( data, "Density_hot", index, NULL );

    index = 0;
    for ( i=0; i<N_Gas; i++ ) {
        if ( !InSlice(i) )
            continue;
        if ( SphP[i].Temp < 1e7 && SphP[i].Temp >= 1e5 ) {
            data[3*index] = P[i].Pos[proj_i] - SliceStart[proj_i];
            data[3*index+1] = P[i].Pos[proj_j] - SliceStart[proj_j];
            data[3*index+2] = SphP[i].Density / Time3 / RhoBaryon;
            index ++;
        }
    }
    field_slice( data, "Density_warm", index, NULL );

    index = 0;
    for ( i=0; i<N_Gas; i++ ) {
        if ( !InSlice(i) )
            continue;
        if ( ( SphP[i].Density / Time3 / RhoBaryon ) < 1e3 && SphP[i].Temp < 1e5 ) {
            data[3*index] = P[i].Pos[proj_i] - SliceStart[proj_i];
            data[3*index+1] = P[i].Pos[proj_j] - SliceStart[proj_j];
            data[3*index+2] = SphP[i].Density / Time3 / RhoBaryon;
            index ++;
        }
    }
    field_slice( data, "Density_cool", index, NULL );

    index = 0;
    for ( i=0; i<N_Gas; i++ ) {
        if ( !InSlice(i) )
            continue;
        if ( ( SphP[i].Density / Time3 / RhoBaryon ) >= 1e3 && SphP[i].Temp < 1e5 ) {
            data[3*index] = P[i].Pos[proj_i] - SliceStart[proj_i];
            data[3*index+1] = P[i].Pos[proj_j] - SliceStart[proj_j];
            data[3*index+2] = SphP[i].Density / Time3 / RhoBaryon;
            index ++;
        }
    }
    field_slice( data, "Density_condensed", index, NULL );

    myfree( data );
#endif
}
