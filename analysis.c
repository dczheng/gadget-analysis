#include "allvars.h"

#define make_group_output_filename( buf, nstr, group_index ) \
    sprintf( buf, "%s/%s/%s_%s_%03i_%.2f_%04i_%c.dat",\
            All.GroupDir, nstr, All.FilePrefix, nstr, ThisTask,\
            All.RedShift, group_index, All.Sproj );

void init_analysis() {
    writelog( "initialize analysis...\n" );
    All.proj_k = All.ProjectDirection;
    All.proj_i = ( All.ProjectDirection + 1 ) % 3;
    All.proj_j = ( All.ProjectDirection + 2 ) % 3;
    All.Sproj = All.ProjectDirection + 'x';
    All.PicSize2 = SQR( All.PicSize );
    slice();
    if ( All.KernelInterpolation )
        init_kernel_matrix();
    init_conv_kernel();
    init_img();
    writelog( "initialize analysis... done.\n" );
    put_block_line;
}

void free_analysis() {
    writelog( "free analysis ...\n" );
    if ( All.KernelInterpolation )
        free_kernel_matrix();
    free_conv_kernel();
    free_img();
    writelog( "free analysis ... done.\n" );
    put_block_line;
}

void gas_density_slice() {
    int num, i;
    char buf[100];
    writelog( "gas density silce ...\n" );
    num = All.SliceEnd[0] - All.SliceStart[0];
    mymalloc2( image.data, sizeof( double ) * num );
    mymalloc2( image.img, sizeof( double ) * SQR( All.PicSize ) );

    for ( i=All.SliceStart[0]; i<num; i++ ) {
        image.data[i] = SphP[i].Density;
    }

    create_dir( "./gas_density" );
    sprintf( buf, "./gas_density/%s_gas_density_%03i_%.2f.dat", All.FilePrefix, ThisTask, All.RedShift );

    make_slice_img( 0 );

    for ( i=0; i<SQR(All.PicSize); i++ ) {
        image.img[i] *= All.UnitMass_in_g / pow( All.UnitLength_in_cm, 2 );
    }


    write_img2( buf, "gas density slice" );
    myfree( image.data );
    myfree( image.img );

    writelog( "gas density silce ... done.\n" );
    put_block_line;
}

int compare_gas_rho( const void *a, const void *b ){
    return ((((struct Sph_Particle_Data* )a)->Density) < (((struct Sph_Particle_Data*)b)->Density)) ? 1: -1;
}

void sort_gas_rho(){
    int i;
    qsort( (void*)SphP, N_Gas, sizeof( struct Sph_Particle_Data ), &compare_gas_rho );
    for ( i=0; i<10; i++ ) {
        printf( "%g\n", SphP[i].Density * ( g2c.g / CUBE(g2c.cm) ) );
    }
    printf( "\n" );
    for ( i=0; i<10; i++ ) {
        printf( "%g\n", SphP[N_Gas-10+i].Density * ( g2c.g / CUBE(g2c.cm) ) );
    }
}

void vel_value() {

    writelog( "velocities value analysis ...\n" );
    FILE *fd;
    char buf[20];
    int i;
    double v;
    sprintf( buf, "vel_%i.txt", ThisTask );
    fd = fopen( buf, "w" );
    for ( i=0; i<N_Gas; i++ ) {
        v = sqrt( pow( P[i].Vel[0], 2 ) + pow( P[i].Vel[1], 2 ) + pow( P[i].Vel[2], 2 ) );
        if ( v > 1000 ) {
            fprintf( fd, "%i %g\n", P[i].ID, v );
        }
    }
    fclose( fd );
    writelog( "velocities value analysis ... done\n" );
    put_block_line;
}

/*
void compute_radio( double v ) {

    mymalloc( img, sizeof( double ) * SQR( PicSize ) );
    memset( img, 0, sizeof( double ) * SQR( PicSize ) );

    dy = dx = BoxSize / PicSize;
    h = All.SofteningTable[0] * ( g2c.cm );
    V = 4.0 / 3.0 * PI * pow( h, 3 );
    ang = h / ang_dis / PI * 180.0 * 60;
    beam = pow( 10.0 / 60, 2.0 ) / ( 4.0 * log(2) );
    img_max = DBL_MIN;
    img_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        xi = (int)( P[i].Pos[0] / dx );
        yi = (int)( P[i].Pos[1] / dy );
        tmp = img[xi * PicSize + yi] += SphP[i].P * V / ( 4*PI*pow(lum_dis,2) ) / SQR(ang)
            * beam * 1e25;
        //tmp = img[xi * PicSize + yi] += SphP[i].P * V / ( 4*PI*pow(lum_dis,2) ) / SQR(ang)
        //    * 1e25;
        img_max = ( tmp > img_max ) ? tmp : img_max;
        img_min = ( tmp < img_min && tmp > 0 ) ? tmp : img_min;
    }
}
*/

void magnetic() {
    long i;
    double B, Bmin, Bmax;

    Bmax = -DBL_MAX;
    Bmin = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        B = sqrt(
                SQR( SphP[i].B[0] ) +
                SQR( SphP[i].B[1] ) +
                SQR( SphP[i].B[2] ) );
        Bmin = vmin( B, Bmin, 0 );
        Bmax = vmax( B, Bmax );
    }
    printf( "Bmin: %g, Bmax: %g\n",
            Bmin, Bmax );

}

void output_rho() {
    long i;
    FILE *fd;
    double rho_max, rho_min;

    fd = fopen( "./rho.txt", "w" );
    rho_max = -DBL_MAX;
    rho_min = DBL_MAX;
    for ( i=0; i<N_Gas; i++ ) {
        fprintf( fd, "%g %g %g %g\n",
                P[i].Pos[0],
                P[i].Pos[1],
                P[i].Pos[2],
                SphP[i].Density );
        rho_max = ( SphP[i].Density > rho_max ) ? SphP[i].Density : rho_max;
        rho_min = ( SphP[i].Density < rho_min ) ? SphP[i].Density : rho_min;
    }
    printf( "rho_max = %g, rho_min = %g \n", rho_max, rho_min );
    fclose( fd );

}

void compute_temperature() {
    double yhelium, u, ne, mu, XH;
    int i;

    writelog( "compute gas temprature...\n" );
    XH = HYDROGEN_MASSFRAC;
    yhelium = ( 1 - XH ) / ( 4 * XH );
    for ( i=0; i<N_Gas; i++ ) {
        u = SphP[i].u * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
        ne = SphP[i].elec;
        mu = ( 1 + 4 * yhelium ) / ( 1 + yhelium + ne );
        SphP[i].Temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
    }
    writelog( "compute gas temprature... done.\n" );
    put_block_line;

}

void gas_state() {

    double TempMin, TempMax, DensMin, DensMax,
           LogTempMin, LogTempMax, LogDensMin, LogDensMax,
           *img, DLogTemp, DLogDens, LogDens, LogTemp, sum,
           GlobLogTempMin, GlobLogTempMax, GlobLogDensMin, GlobLogDensMax;
    int i, j, k, PicSize;
    char buf[200];

    TempMin = DensMin = DBL_MAX;
    TempMax = DensMax = DBL_MIN;
    writelog( "plot gas state...\n" );

    //PicSize_tmp = All.PicSize;
    PicSize = All.PicSize;


    mymalloc2( img, sizeof( double ) * SQR( PicSize ) );

    for ( i=0; i<N_Gas; i++ ) {
        TempMin = vmin( SphP[i].Temp, TempMin, 1 );
        TempMax = vmax( SphP[i].Temp, TempMax );
        DensMin = vmin( SphP[i].Density/All.RhoBaryon, DensMin, 1 );
        DensMax = vmax( SphP[i].Density/All.RhoBaryon, DensMax );
    }

    LogDensMin = log10( DensMin );
    LogDensMax = log10( DensMax );
    LogTempMin = log10( TempMin );
    LogTempMax = log10( TempMax );
    find_global_value( LogDensMin, GlobLogDensMin, MPI_DOUBLE, MPI_MIN );
    find_global_value( LogDensMax, GlobLogDensMax, MPI_DOUBLE, MPI_MAX );
    find_global_value( LogTempMin, GlobLogTempMin, MPI_DOUBLE, MPI_MIN );
    find_global_value( LogTempMax, GlobLogTempMax, MPI_DOUBLE, MPI_MAX );

    //LogTempMax = 7;

    DLogTemp = ( LogTempMax - LogTempMin ) / PicSize;
    DLogDens = ( LogDensMax - LogDensMin ) / PicSize;
    writelog( "DensMin: %g, DensMax: %g, TempMin: %g, TempMax: %g\n"
            "LogDensMin: %g, LogDensMax: %g, LogTempMIn: %g, LogTempMax: %g\n"
            "GlobLogDensMin: %g, GlobLogDensMax: %g,\nGlobLogTempMIn: %g, GlobLogTempMax: %g\n",
            DensMin, DensMax, TempMin, TempMax,
            LogDensMin, LogDensMax, LogTempMin, LogTempMax,
            GlobLogDensMin, GlobLogDensMax, GlobLogTempMin, GlobLogTempMax );

    ZSPRINTF( 0, "1" );
    for ( k=0; k<N_Gas; k++ ) {
        LogTemp = SphP[k].Temp;
        if ( LogTemp < 0 )
            endrun();
        LogTemp = ( LogTemp == 0 ) ? LogTempMin : log10( LogTemp );

        //if ( LogTemp > LogTempMax )
        //    continue;

        LogDens = SphP[k].Density / All.RhoBaryon;
        if ( LogDens < 0 )
            endrun();
        LogDens = ( LogDens == 0 ) ? LogDensMin : log10( LogDens );

        i = ( LogTemp - LogTempMin ) / DLogTemp;
        check_picture_index( i );

        j = ( LogDens - LogDensMin ) / DLogDens;
        check_picture_index( j );
        if ( i < 0 || i >= All.PicSize || j < 0 || j >= All.PicSize)
            printf( "%i, %i\n", i, j );
        img[ i*PicSize + j ]++;
     //  printf( "%g ", img[ i*PicSize +j ] ) ;
    }
    ZSPRINTF( 0, "2" );

    for ( i=0, sum=0; i<SQR(PicSize); i++ )
        sum += img[i];

    for ( i=0; i<SQR(PicSize); i++ )
        img[i] /= sum;

    create_dir( "./gas_state" );
    sprintf( buf, "./gas_state/%s_gas_state_%03i_%.2f.dat", All.FilePrefix, ThisTask, All.RedShift );

    ZSPRINTF( 0, "3" );
    image.img = img;
    img_xmin = LogDensMin;
    img_xmax = LogDensMax;
    img_ymin = LogTempMin;
    img_ymax = LogTempMax;
    write_img2( buf, "gas state" );
    myfree( img );

    //All.PicSize = PicSize_tmp;
    put_block_line;

}

void gas_temperature_slice() {

    int num, i, PicSize2;
    char buf[100];
    writelog( "gas temperature silce ...\n" );
    PicSize2 = All.PicSize2;
    num = All.SliceEnd[0] - All.SliceStart[0];
    mymalloc2( image.data, sizeof( double ) * num );
    mymalloc2( image.img, sizeof( double ) * SQR( All.PicSize ) );

    for ( i=All.SliceStart[0]; i<num; i++ ) {
        image.data[i] = SphP[i].Temp;
    }

    create_dir( "./gas_temperature" );
    sprintf( buf, "./gas_temperature/%s_gas_temperature_%03i_%.2f.dat", All.FilePrefix, ThisTask, All.RedShift );

    make_slice_img( 0 );

    for ( i=0; i<PicSize2; i++ ){
        image.img[i] *= All.UnitLength_in_cm;
    }

    write_img2( buf, "gas temperature slice" );

    myfree( image.data );
    myfree( image.img );
    writelog( "gas Temperature silce ... done.\n" );
    put_block_line;

}

void mass_function() {

    int g, *num, i, i_m_min, i_m_max;
    double m_min, m_max;
    char buf[100];
    FILE *fd;
    writelog( "compute mass function ...\n" );

    m_min = DBL_MAX;
    m_max = -DBL_MAX;

    for ( g=0; g<Ngroups; g++ ) {
        m_min = vmin( m_min, Gprops[g].mass, 0 );
        m_max = vmax( m_max, Gprops[g].mass );
    }

    m_min *= 1e10;
    m_max *= 1e10;

    writelog( "m_min: %g, m_max: %g\n", m_min, m_max );

    i_m_min = floor( log10(m_min) );
    i_m_max = floor( log10(m_max) );

    writelog( "i_m_min: %i, i_m_max: %i\n", i_m_min, i_m_max );

    if ( i_m_min == i_m_max ) {
        endrun( "some error appear!" );
    }

    mymalloc2( num, sizeof(int) * ( i_m_max - i_m_min + 1 ) );

    for ( g=0; g<Ngroups; g++ ) {
        num[ ( int ) ( floor( log10( Gprops[g].mass * 1e10 ) ) - i_m_min ) ] ++;
    }

    for ( i=i_m_max-1; i>=i_m_min; i-- ) {
        num[i-i_m_min] += num[i-i_m_min+1];
    }

    /*
    for ( i=m_min; i<=m_max; i++ ) {
        writelog( "%i\n", num[i] );
    }
    */

    create_dir( "./mf" );
    sprintf( buf, "./mf/%s_mf_%.2f.dat", All.FilePrefix, All.RedShift );

    fd = fopen( buf, "w" );

    for ( i=i_m_min; i<=i_m_max; i++ ) {

        fprintf( fd, "%i %g\n", i, num[i-i_m_min] / CUBE( All.BoxSize / All.MpcFlag ) );

    }

    fclose( fd );

    myfree( num );

    writelog( "compute mass function ... done.\n" );
    put_block_line;
}

int group_present( long index ) {

    if ( ( index >= All.GroupIndexMin ) &&
         ( index <= All.GroupIndexMax ) &&
         ( Gprops[index].mass * 1e10 >= All.GroupMassMin) )
            return 1;
    else {

        return 0;
    }

}

double particle_radio( double v, long i ) {

#define BCMB0 (3.24e-6) // Gauss
    double C, B, Ub, vl, P;

    C = SphP[i].CRE_C;
    //printf( "%g\n", C );
    C = C * SphP[i].Density / ( ELECTRON_MASS /(g2c.g) );

    C /= CUBE( g2c.cm );

    B = sqrt( pow( SphP[i].B[0], 2 ) + pow( SphP[i].B[1], 2 ) + pow( SphP[i].B[2], 2 ) );
    Ub = B * B / ( 8 * PI );
    Ub += SQR( BCMB0 ) * pow( All.Time, -4 ) / ( 8*PI );

    vl = ELECTRON_CHARGE * B / ( 2 * PI * ELECTRON_MASS * LIGHT_SPEED );

    P = C * 2.0 / 3.0 * LIGHT_SPEED * Ub * THOMSON_CROSS_SECTION *
        pow( v / vl-1, (1-SphP[i].CRE_Alpha) / 2 ) / vl;
    //printf( "%g\n", P );
    //
    P = P * ( 4.0/3.0 * PI * CUBE( All.SofteningTable[0] * g2c.cm ) );

    // Unit: erg / Hz / s

    return P;

}

double group_luminosity( double nu, long index ) {

    long p;
    double F;
    struct group_properties *g;

    g = &Gprops[index];

    p = g->Head;
    F = 0;

    while( p >= 0 ) {

        if ( P[p].Type == 0 )
            F += particle_radio( nu, p );
        p = FoFNext[p];

    }

    return F;

}

void group_flux( double nu, long index, double *flux, double *flux_nosr ) {

    double com_dis, lum_dis, L;
    struct group_properties *g;

    g = &Gprops[index];
    com_dis = comoving_distance( All.Time );
    lum_dis = luminosity_distance( All.Time ) * g2c.cm;

    L = group_luminosity( nu, index );
    *flux_nosr = L / ( 4.0 * PI * SQR( lum_dis ) );
    *flux = *flux_nosr / ( SQR( g->size * 2  / com_dis ) );

}

void group_spectrum() {

    long index;
    int vN, i, signal;
    double *v, *flux, vmin, vmax, dv, *flux_nosr;
    char buf[100];
    FILE *fd1, *fd2;

    writelog( "group spectrum ...\n" );

    vN = All.NuNum;
    signal = vN / 10;
    vmin = All.NuMin;
    vmax = All.NuMax;
    dv = log( vmax/vmin) / vN;

    mymalloc1( v, sizeof( double ) * vN );
    mymalloc1( flux, sizeof( double ) * vN );
    mymalloc1( flux_nosr, sizeof( double ) * vN );

    sprintf( buf, "%s/Spectrum", All.GroupDir );
    create_dir( buf );

    sprintf( buf, "%s/Spectrum/%s_%03i_spec_%.2f.dat",
            All.GroupDir, All.FilePrefix, ThisTask, All.RedShift );
    fd1 = fopen( buf, "w" );

    sprintf( buf, "%s/Spectrum/%s_%03i_spec_%.2f_nosr.dat",
            All.GroupDir, All.FilePrefix, ThisTask, All.RedShift );
    fd2 = fopen( buf, "w" );

    fprintf( fd1, "0  0  " );
    fprintf( fd2, "0  0  " );

    for ( i=0; i<vN; i++ ) {
        v[i] = exp(log(vmin) + i * dv);
        fprintf( fd1, "%g  ", v[i] );
        fprintf( fd2, "%g  ", v[i] );
    }
    fprintf( fd1, "\n" );
    fprintf( fd2, "\n" );

    for ( index=0; index<Ngroups; index++ ) {

        if ( !group_present( index ) )
            break;

        for ( i=0; i<vN; i++ ) {

            if ( ( i % signal == 0 ) || ( i == vN - 1 ) )
                writelog( "group: %li, [%i]: %g MHz ...\n", index, i, v[i] );

            group_flux( v[i] * 1e6, index, flux+i, flux_nosr+i );
              //p[i] = group_luminosity( v[i] * 1e6, index );
        }



        fprintf( fd1, "%li  %g  ", index, Gprops[index].mass * 1e10 );
        fprintf( fd2, "%li  %g  ", index, Gprops[index].mass * 1e10 );

        for ( i=0; i<vN; i++ ) {
            fprintf( fd1, "%g  ", flux[i] );
            fprintf( fd2, "%g  ", flux_nosr[i] );
        }
        fprintf( fd1, "\n" );
        fprintf( fd2, "\n" );

    }

    fclose( fd1 );
    fclose( fd2 );

    myfree( v );
    myfree( flux );
    myfree( flux_nosr );

    writelog( "group spectrum ... done.\n" );

}

void group_spectrum_index() {

    double *spec, *v, vmin, vmax, dv, cov00, cov01, cov11, c0,
           *spec_index, *spec_index_err, nu, L, dL, *mass;
    int vN, k, PicS, PicS2, ii, jj, i, j, x, y, xo, yo, flag,
        index, p;
    struct group_properties g;
    char buf[100],
         *spec_index_str="Spectrum_index",
         *spec_index_err_str="Spectrum_index_err";

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

        L = g.size;
        dL = 2 * L / PicS;

        for ( k=0; k<vN; k++ ) {

            nu = pow( 10, v[k] ) * 1e6;
            p = g.Head;
            for ( i=0; i<g.Len; i++, p=FoFNext[p] ) {

                if ( P[p].Type != 0 )
                    continue;

                ii = PERIODIC( P[p].Pos[x] - g.cm[x] ) / dL + xo;
                jj = PERIODIC( P[p].Pos[y] - g.cm[y] ) / dL + yo;
                check_picture_index( ii );
                check_picture_index( jj );
                spec[ii*PicS*vN + jj * vN + k] += particle_radio( nu, p );

            }
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
        sprintf( buf, "%s/%s/%s_%s_%03i_%.2f_%04i_%c.dat",
            All.GroupDir, spec_index_str, All.FilePrefix, spec_index_err_str, ThisTask,
            All.RedShift, index, All.Sproj );
        write_img1( buf, spec_index_err_str );

    }

    myfree( v );
    myfree( spec );
    myfree( spec_index );
    myfree( spec_index_err );
    writelog( "group spectrum index... done.\n" );

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
        case GROUP_HGEN:
            if ( All.GroupHgen )
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
        case GROUP_HGEN:
            strcpy( buf, "Hge_n" );
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
            case GROUP_HGEN:
                if ( All.ReadHge == 0 ) {
                    printf( "`ReadHge` is required by `GroupHgen`.\n" );
                    flag = 1;
                }
                break;
            case GROUP_RAD:
                if ( All.ReadHge == 0  || All.ReadB == 0) {
                    printf( "`ReadHge` and `RaadB` is required by `GroupRad`.\n" );
                    flag = 1;
                }
                break;
        }
    }

    if ( flag )
        endrun( 20180926 );

}


void group_analysis() {

    long *npart;
    struct group_properties g;
    double L, dL, *mass, B, PP, nu, n,
            lum_dis, com_dis;// com_dis, lum_dis, h, ang;
    int *num, g_index, i, j, p, PicSize, x, y,
         xo, yo, PicSize2, pic_index, ii, jj;

    char buf[100], buf1[100];
    double *data[GROUP_FILED_NBLOCKS];

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
        L = g.size;

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
            ii = PERIODIC( P[p].Pos[x] - g.cm[x] ) / dL + xo;
            jj = PERIODIC( P[p].Pos[y] - g.cm[y] ) / dL + yo;
            //printf( "%li, %li\n", ii, jj );
            check_picture_index( ii );
            check_picture_index( jj );
            pic_index= ii*PicSize + jj;

            num[ pic_index ] ++;

            if ( group_filed_present( GROUP_DENS ) )
                data[GROUP_DENS][pic_index]       += SphP[p].Density;

            if ( group_filed_present( GROUP_TEMP ) )
                data[GROUP_TEMP][pic_index]       += SphP[p].Temp * SphP[p].Density;

            if ( group_filed_present( GROUP_SFR ) )
                data[GROUP_SFR][pic_index]        += SphP[p].sfr * SphP[p].Density;

            if ( group_filed_present( GROUP_MAG ) ) {
                B = sqrt( SQR( SphP[p].B[0] ) +
                          SQR( SphP[p].B[1] ) +
                          SQR( SphP[p].B[2] ) );
                B *= 1e6;  // convert G to muG
                data[GROUP_MAG][pic_index]        += B * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_MACH ) )
                data[GROUP_MACH][pic_index]       += SphP[p].MachNumber * SphP[p].Density;

            if ( group_filed_present( GROUP_HGEN ) ) {
                n = SphP[p].CRE_n *  SphP[p].Density / ( ( ELECTRON_MASS / ( g2c.g ) ) );
                n /= CUBE( g2c.cm );
                data[GROUP_HGEN][pic_index]      += n * SphP[p].Density;
            }

            if ( group_filed_present( GROUP_HGEN ) ) {
                nu = All.Freq * 1e6;
                PP = particle_radio( nu, p );
                data[GROUP_RAD][pic_index]        += PP * SphP[p].Density;
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
            data[GROUP_DENS][i]  *= (( g2c.g ) / CUBE( g2c.cm ) / CUBE( All.Time ));

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

    put_block_line;

    if ( All.GroupSpec )
        group_spectrum();

    put_block_line;

    if ( All.GroupSpecIndex )
        group_spectrum_index();

    put_block_line;

}

void total_radio_spectrum() {

    long index;
    int Nnu, i, signal, j, k;
    double *nu, *flux, numin, numax, dnu, tmp,
           *alist, *fluxlist, *dislist;
    char buf[100];
    FILE *fd;

    writelog( "total radio spectrum ...\n" );

    Nnu = All.NuNum;
    signal = Nnu / 10;
    numin = All.NuMin;
    numax = All.NuMax;
    dnu = log( numax/numin) / Nnu;

    mymalloc1( nu, sizeof( double ) * Nnu );
    mymalloc2( flux, sizeof( double ) * Nnu );

    for ( i=0; i<Nnu; i++ ) {
        nu[i] = exp(log(numin) + i * dnu);
    }

    if ( All.RedShift > 1e-10 )
        for ( i=0; i<Nnu; i++ ) {
            if ( ( i % signal == 0 ) || ( i == Nnu - 1 ) )
                writelog( "total spectrum [%i]: %g MHz ...\n", i, nu[i] );

            for ( index=0; index<N_Gas; index++ ) {
                flux[i] += particle_radio( nu[i]*1e6 / All.Time, index );
            }
        }

    tmp = 1.0 / ( All.Time * ( 4.0 * PI * SQR( All.ComDis * g2c.cm ) ) * ( SQR( All.BoxSize / All.ComDis ) ) );

    for ( i=0; i<Nnu; i++ )
        flux[i] *= tmp;

    create_dir( "TotalSpec" );

    MPI_Barrier( MPI_COMM_WORLD );
    mymalloc1( alist, sizeof( double ) * NTask );
    mymalloc1( dislist, sizeof( double ) * NTask );
    mymalloc1( fluxlist, sizeof( double ) * NTask * Nnu );

    MPI_Gather( &All.Time, 1, MPI_DOUBLE, alist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Gather( flux, Nnu, MPI_DOUBLE, fluxlist, Nnu, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    if ( ThisTask == 0 ) {
        for ( i=0; i<NTask-1; i++ )
            for ( j=i; j<NTask; j++ ) {
                if ( alist[i] < alist[j] ) {

                        tmp = alist[i];
                        alist[i] = alist[j];
                        alist[j] = tmp;

                        for ( k=0; k<Nnu; k++ ) {
                            tmp = fluxlist[ i*Nnu + k ];
                            fluxlist[ i*Nnu + k ] = fluxlist[ j*Nnu + k ];
                            fluxlist[ j*Nnu + k ] = tmp;
                        }

                }
            }

        for ( i=0; i<NTask-1; i++ )
            dislist[i] = comoving_distance( alist[i] );

        for ( i=0; i<NTask*Nnu; i++ )
            fluxlist[i] /= All.BoxSize;

        memset( flux, 0, sizeof(double)*Nnu );

        for ( i=0; i<Nnu; i++ )
            for ( j=0; j<NTask-1; j++) {
                flux[i] += 0.5 * ( fluxlist[j*Nnu+i] + fluxlist[(j+1)*Nnu+i] ) * ( dislist[j+1] - dislist[j] );
            }

        sprintf( buf, "./TotalSpec/%s_Spec_Tot.dat", All.FilePrefix );
        fd = fopen( buf, "w" );

        for ( i=0; i<Nnu; i++ )
            fprintf( fd, "%g %g\n", nu[i], flux[i] );
        fclose( fd );
    }

    myfree( alist );
    myfree( dislist );
    myfree( fluxlist );
    myfree( nu );
    myfree( flux );

    writelog( "total radio spectrum ... done.\n" );

}

void analysis(){
    init_analysis();

    if ( (All.GasTemperature ||
          All.GasState ||
          All.Group )
         && All.ReadTemp == 0  ) {
        compute_temperature();
    }

    if ( All.GasState ) {
        gas_state();
    }

    if ( All.GasDensity )
        gas_density_slice();

    if ( All.GasTemperature )
        gas_temperature_slice();

    if ( All.FoF ) {
        if ( All.FoFRead )
            fof_read();
        else
            fof();
    }

    if ( All.MF )
        mass_function();

    if ( All.Group ) {

        if ( All.FoF == 0 )
            endrun( "FoF is required by Group Analysis!" );

        create_dir( All.GroupDir );
        group_analysis();

    }

    if ( All.TotSpec ) {
        total_radio_spectrum();

    if( All.FoF )
        fof_free();
    }

    //tree_build();
    //tree_free();
    //output_rho();
    //vel_value();
    //sort_gas_rho();
    free_analysis();
}
