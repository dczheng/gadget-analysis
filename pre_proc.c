#include "allvars.h"

struct sphp_sort{
    ParticleData *P;
    SphParticleData *SphP;
};

void find_id() {
    int bits;
    long i, num, offset;
    mytimer_start();
    writelog( "find id ...\n" );

    for ( bits=0; GENERATIONS > (1<<bits); bits++ );
    writelog( "bits: %i\n", bits );
    writelog( "find gas id ...\n" );
    for ( i=0; i<N_Gas; i++ ){
        P[i].ID <<= bits;
        P[i].ID >>= bits;
    }

    num = get_particle_num( 4 );
    offset = get_particle_offset( 4 );
    if ( num != 0 ) {
        writelog( "Star Particle Offset: %li\n", offset );
        writelog( "find star id ...\n" );
        for ( i=0; i<num; i++ ) {
            P[offset+i].ID <<= bits;
            P[offset+i].ID >>= bits;
        }
    }

    num = get_particle_num( 5 );
    offset = get_particle_offset( 5 );
    if ( num != 0 ) {
        writelog( "Black Hole Particle Offset: %li\n", offset );
        writelog( "find Black Hole id ...\n" );
        for ( i=0; i<num; i++ ) {
            P[offset+i].ID <<= bits;
            P[offset+i].ID >>= bits;
        }
    }

    writelog( "find id ... done.\n" );
    mytimer_end();
    put_block_line;
}

void construct_id_to_index() {
    long idmax, idmin, i, idn, N;

    mytimer_start();
    writelog( "construct id to index ...\n" );
    idmax = -1;
    idmin = LONG_MAX;

    for ( i=0; i<NumPart; i++ ) {
        idmax = vmax( idmax, P[i].ID );
        idmin = vmin( idmin, P[i].ID, 0 );
    }

    idn = idmax - idmin + 1;
    writelog( "Min ID: %li, Max ID: %li, ID Num: %li\n", idmin, idmax, idn );

    mymalloc3( id_to_index, idn * sizeof( long ), -1 );

    id_to_index--;

    N = N_Gas + get_particle_num( 1 );

    for ( i=0; i<N; i++ ) {
        id_to_index[ P[i].ID ] = i;
    }

    /*
    long offset, num, id;
    offset = get_particle_offset( 4 );
    num = get_particle_num( 4 );
    for ( i=offset; i<offset+num; i++ ) {
        id = P[i].ID;
        if ( P[id_to_index[id]].ID != id )
            printf( "%li %li %li\n", id, id_to_index[id], P[id_to_index[id]].ID );
    }
    */

    writelog( "construct id to index ... done.\n" );
    mytimer_end();

    put_block_line;
}

int compare_pos( const void *a, const void *b ) {

    double *ap, *bp;

    ap = ( ( ParticleData* )a ) -> Pos;
    bp = ( ( ParticleData* )b ) -> Pos;

    if ( ap[0] > bp[0] )
        return 1;

    if ( ap[0] < bp[0] )
        return -1;

    if ( ap[1] > bp[1] )
        return 1;

    if ( ap[1] < bp[1] )
        return -1;

    if ( ap[2] > bp[2] )
        return 1;

    if ( ap[2] < bp[2] )
        return -1;

    return 0;

}

int compare_sphp_pos( const void *a, const void *b ) {

    ParticleData *aa, *bb;

    aa = ( ( struct sphp_sort* ) a ) -> P;
    bb = ( ( struct sphp_sort* ) b ) -> P;

    return compare_pos( aa, bb );

}

void sort_particle_by_pos() {

    long i, n, num, offset;

    char *flag;

    ParticleData p;
    SphParticleData sp;

    struct sphp_sort *ss;


    writelog( "sort particle ...\n" );
    mytimer_start();

    writelog( "sort sph particle ...\n" );

    mymalloc1( ss, sizeof( struct sphp_sort ) * N_Gas );
    mymalloc2( flag, N_Gas );

    for( i=0; i<N_Gas; i++ ) {

        ss[i].P = &P[i];
        ss[i].SphP = &SphP[i];

    }

    qsort( ss, N_Gas, sizeof( struct sphp_sort ), &compare_sphp_pos );

    for( i=0; i<N_Gas; i++ ) {

        if ( flag[i] == 0 ) {

            p = P[i];
            sp = SphP[i];

            n = i;

            while( (ss[n].P - P) != i ) {

                flag[n] = 1;

                P[n] = *ss[n].P;
                SphP[n] = *ss[n].SphP;

                n = ss[n].P - P;

            }

            flag[n] = 1;

            P[n] = p;
            SphP[n] = sp;

        }

    }

    /*
    put_block_line;

    for( i=0; i<100; i++ )  {
        printf( "[%li] %e %e %e\n",
                i,
                P[i].Pos[0],
                P[i].Pos[1],
                P[i].Pos[2] );
    }

    for( i=0; i<100; i++ )  {
        printf( "[%li] %e %e %e\n",
                N_Gas - 100 + i,
                P[N_Gas - 100 + i].Pos[0],
                P[N_Gas - 100 + i].Pos[1],
                P[N_Gas - 100 + i].Pos[2] );
    }

    put_block_line;
    */

    myfree( flag );
    myfree( ss );

    writelog( "sort sph particle ... done.\n" );
    mytimer();

    for( i=1; i<6; i++ ) {

        offset = get_particle_offset( i );
        num = get_particle_num( i );

        if ( num>0 ) {
            writelog( "sort `%li` particle ... \n", i );
            qsort( P+offset, num, sizeof( ParticleData ), &compare_pos );
            writelog( "sort `%li` particle ... done.\n", i );
            mytimer();
        }

    }

    /*
    for( i=0; i<NumPart-1; i++ ) {
        if ( compare_pos( &P[i], &P[i+1]) &&
                P[i].Type == P[i+1].Type ) {
            printf( "[%li] ( %e, %e, %e )\n"
                    "[%li] ( %e, %e, %e )\n",
                    i, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
                    i+1, P[i+1].Pos[0], P[i].Pos[1], P[i].Pos[2] );

            endruns( "Error 20181109" );
        }
    }
    */

    mytimer_end();
    writelog( "sort particle ...done.\n" );

    put_block_line;

}

void test_sort() {

    ParticleData *Ptmp;
    SphParticleData *SphPtmp;

    long i, li, ri, mi;
    int signal;

    writelog( "test sort ...\n" );
    mymalloc1( Ptmp, sizeof( ParticleData ) * N_Gas );
    mymalloc1( SphPtmp, sizeof( SphParticleData ) * N_Gas );

    memcpy( Ptmp, P, sizeof( ParticleData ) * N_Gas );
    memcpy( SphPtmp, SphP, sizeof( SphParticleData ) * N_Gas );

    sort_particle_by_pos();

    mytimer_start();

    signal = N_Gas / 100;

    for( i=0; i<N_Gas; i++ ) {

        if ( i % signal == 0 ) {
            writelog( "[%6.2f%%] [%li] [%li]\n", (double)i / N_Gas * 100, N_Gas, i );
        }

        li = 0;
        ri = N_Gas-1;
        mi = ( li+ri ) / 2;

        /*
        if ( (double)i / N_Gas > 0.68 &&
                (double)i / N_Gas < 0.70 )
            writelog( "%li\n", i );
            */

        if ( compare_pos( &Ptmp[i], &P[ri] ) == 0 ||
                compare_pos( &Ptmp[i], &P[li] ) == 0 ) {

            if ( compare_pos( &Ptmp[i], &P[ri] ) == 0 )
                mi = ri;
            else
                mi = li;
        }

        else {

            while( compare_pos( &Ptmp[i], &P[mi] ) != 0 ) {

                if ( li == ri )
                    endruns( "Error 20181108 0" );

                if ( compare_pos( &Ptmp[i], &P[mi] ) > 0 )
                    li = mi;
                else
                    ri = mi;

                if ( li > ri )
                    endruns( "Error 20181108 1" );

                mi = ( li+ri ) / 2;

            }

        }

        /*
        if ( Ptmp[i].ID != P[mi].ID ||
                SphPtmp[i].Density != SphP[mi].Density )
                */
        if ( memcmp( &Ptmp[i], &P[mi], sizeof( ParticleData ) ) != 0 ||
                memcmp( &SphPtmp[i], &SphP[mi], sizeof( SphParticleData ) ) != 0 )
            endruns( "Error 20181108 2" );

    }

    myfree( Ptmp );
    myfree( SphPtmp );
    writelog( "sort is ok ...\n" );
    mytimer_end();

    writelog( "test sort ... done.\n" );

}

void merge_particle( int pt ) {

    long offset, num0, num1, i, i1, i2, j, n;

    if ( pt == 0 ) {
        writelog( "Error: particle `0` can not be merge !!!\n" );
        endrun( 20181109 );
    }

    writelog( "merge particle `%i` ...\n", pt );

    offset = get_particle_offset( pt );
    num0 = get_particle_num( pt );

    writelog( "offset: %li, num: %li\n", offset, num0 );

    num1 = num0;

    for( i=offset; i<offset+num1-1; i++ ) {

        if ( compare_pos( &P[i], &P[i+1] ) == 0 ) {

            i1 = i2 = i;
            while( compare_pos( &P[i2], &P[i2+1] ) == 0 && i2 < offset+num1-1 )
                i2++;

            for( j=i1+1; j<=i2; j++ ) {
                P[i1].Mass += P[j].Mass;
                P[j] = P[ offset+num1-1 ];
                num1 --;
            }

            //writelog( "[%li]\n", i );
        }
    }

    n = num0 - num1;

    writelog( "merge num: %li\n", n );

    //return;

    if ( n != 0 ) {

        header.npartTotal[pt] = ( ( num1<<32 ) >> 32 );
        header.npartTotalHighWord[pt] = ( num1>>32 );

        for( i=offset+num0; i<NumPart; i++ ) {
            P[offset+num1 + ( i-(offset+num0) )] = P[i];
        }

        NumPart -= n;
    }

    writelog( "merge particle `%i` ... done.\n", pt );
    put_block_line;

}

void attach_particle_to_gas( int pt ) {

    long i, li, ri, mi, offset, num, attach_num;

    if ( (pt != 4) && (pt != 5) ) {
        writelog( "Error: Attach is only for particle `4` or `5`\n" )
        endrun( 20181109 );
    }

    offset = get_particle_offset( pt );
    num = get_particle_num( pt );

    pt -= 4;

    attach_num = 0;

    for ( i=offset; i<offset+num; i++ ) {

        li = 0;
        ri = N_Gas - 1;
        mi = ( li+ri ) / 2;

        if ( compare_pos( &P[i], &P[ri] ) == 0 ||
                compare_pos( &P[i], &P[li] ) == 0 ) {

            if ( compare_pos( &P[i], &P[ri] ) == 0 )
                mi = ri;
            else
                mi = li;
        }

        else {

            while( compare_pos( &P[i], &P[mi] ) != 0 ) {

                if ( compare_pos( &P[i], &P[mi] ) > 0 )
                    li = mi;
                else
                    ri = mi;

                if ( ri - li == 1 )
                    break;

                mi = ( li+ri ) / 2;
            }
        }

        if ( compare_pos( &P[i], &P[mi] ) == 0 ) {

            attach_num ++;

            if ( SphP[mi].Star_BH_Num[pt] == 0 ) {
                SphP[mi].Star_BH_MaxNum[pt] = 4;
                SphP[mi].Star_BH_Index[pt] = malloc( sizeof( MyIDType ) * 4 );
            }

            if ( SphP[mi].Star_BH_Num[pt] == SphP[mi].Star_BH_MaxNum[pt] ) {
                SphP[mi].Star_BH_MaxNum[pt] += 4;
                SphP[mi].Star_BH_Index[pt] = realloc( SphP[mi].Star_BH_Index[pt],
                        sizeof( MyIDType ) * SphP[mi].Star_BH_MaxNum[pt] );
            }

            SphP[mi].Star_BH_Index[pt][ SphP[mi].Star_BH_Num[pt] ] = i;
            SphP[mi].Star_BH_Num[pt] ++;

        }
    }

    if ( attach_num != 0 ) {
        writelog( "Attatch %li `%i` particle to gas.\n", attach_num, pt+4 );
        put_block_line;
    }

}

void test_attach_particle_to_gas() {

    long offset45[2], i;
    int j, ii, n;

    for( i=0; i<2; i++ )
        offset45[i] = get_particle_offset( i+4 );

    ii = 10;
    n = 10;

    for( i=0; i<n; i++ )
        for( j=0; j<3; j++ )
            P[i].Pos[j] = P[offset45[0] + i].Pos[j];

    for( i=ii; i<ii+n; i++ )
        for( j=0; j<3; j++ )
            P[i].Pos[j] = P[offset45[1] + i].Pos[j];

}

void pre_proc() {

    writelog( "pre proc ...\n" )

    find_id();

    //test_sort();

    //test_attach_particle_to_gas();

    sort_particle_by_pos();

    merge_particle( 4 );
    merge_particle( 5 );

    attach_particle_to_gas( 4 );
    attach_particle_to_gas( 5 );

    //endrun( 20181108 );

    writelog( "pre proc ... done.\n" )

    put_block_line;

}

