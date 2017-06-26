#include "allvars.h"

#define BINS 128
#define BINS2 ( 128 * 128 )
#define BINS3 ( 128 * 128 * 128 )

struct io_header header1;
struct Particle_Struct Particle1[2], Particle2[2];
struct Particle_Line_Struct{
    long *index;
    long num;
    long line_len;
};
struct Particle_Line_Struct Particle_Line[2][ BINS3 ];

void Add_Flag_IC_Info( char *fn ) {
    long tmp=0;
    hdf5_file = H5Fopen( fn, H5F_ACC_RDWR, H5P_DEFAULT );
    hdf5_group = H5Gopen( hdf5_file, "/Header" );
    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_IC_Info", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &tmp );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );
    H5Gclose( hdf5_group );
    H5Fclose( hdf5_file );
}

void Disable_PartType3( char *fn ) {
    int npart[6], i;
    hdf5_file = H5Fopen( fn , H5F_ACC_RDWR, H5P_DEFAULT );
    hdf5_group = H5Gopen( hdf5_file, "/Header" );

    hdf5_attribute = H5Aopen( hdf5_group, "NumPart_ThisFile", H5P_DEFAULT );
    H5Aread( hdf5_attribute, H5T_NATIVE_INT, npart );
    for ( i=0; i<6; i++ ) {
        fprintf( stdout, "%i   ", npart[i] );
    }
    npart[3] = 0;
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, npart );
    H5Aclose( hdf5_attribute );

    hdf5_attribute = H5Aopen( hdf5_group, "NumPart_Total", H5P_DEFAULT );
    H5Aread( hdf5_attribute, H5T_NATIVE_INT, npart );
    for ( i=0; i<6; i++ ) {
        fprintf( stdout, "%i   ", npart[i] );
    }
    npart[3] = 0;
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, npart );
    H5Aclose( hdf5_attribute );

    hdf5_attribute = H5Aopen( hdf5_group, "NumPart_Total_HighWord", H5P_DEFAULT );
    H5Aread( hdf5_attribute, H5T_NATIVE_INT, npart );
    for ( i=0; i<6; i++ ) {
        fprintf( stdout, "%i   ", npart[i] );
    }
    npart[3] = 0;
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, npart );
    H5Aclose( hdf5_attribute );

    H5Gclose( hdf5_group );
    H5Fclose( hdf5_file );
}

void Construct_Particle_Line() {
    long i,index,pt;
    float dx, dy, dz, x, y, z;
    fputs( sep_str, stdout );
    fputs( "Construct Particle_Line ...\n", stdout );
    dx = header.BoxSize / BINS;
    dy = header.BoxSize / BINS;
    dz = header.BoxSize / BINS;
    struct Particle_Line_Struct *p;
    for ( pt=0; pt<2; pt++ ) {
        for ( i=0; i<BINS3; i++ )
            Particle_Line[pt][i].line_len = 0;
    }
    for ( pt=0; pt<2; pt++ ) {
        //fprintf( stdout, "%li\n", Particle[pt].num );
        for ( i=0; i<Particle[pt].num; i++ ) {
            x = Particle[pt].pos[i*3+0];
            y = Particle[pt].pos[i*3+1];
            z = Particle[pt].pos[i*3+2];
            index = (int)( x/dx ) * BINS2 + (int)( y/dy ) * BINS + (int)( z/dz );
            p = &( Particle_Line[pt][index] );
            /*
            if ( 1 == pt  && index == 0) {
                fprintf( stdout, "i=%i, x=%f y=%f z=%f, dx=%f dy=%f dz=%f\n",
                        i, x, y, z, dx, dy, dz );
                fprintf( stdout, "num = %li\n", p->num );
            }
            */
            if ( p->line_len == 0 ) {
                p->line_len = 128;
                p->index = ( long* ) malloc( p->line_len * sizeof( long ) );
                p->num = 0;
                p->index[p->num] = i;
                p->num ++;
                continue;
            }
            if ( p->num == p->line_len ) {
                p->line_len += 128;
                p->index = ( long* ) realloc( p->index,
                                    p->line_len * sizeof( long) );
            }
            p->index[p->num] = i;
            p->num ++;
        }
    }
    fputs( sep_str, stdout );
}

void Generate_Mass() {
    long pt, i, j, index;
    float m;
    struct Particle_Line_Struct *p;
    fputs( sep_str, stdout );
    fprintf( stdout, "Generate Mass ...\n", stdout );
    for ( pt=0; pt<2; pt++ ) {
        Particle1[pt].m = ( float* ) malloc( sizeof( float ) * BINS3 );
        for ( i=0; i<BINS3; i++ ) {
            p = &( Particle_Line[pt][i] );
            Particle1[pt].m[i] = 0;
            for ( j=0; j<p->num; j++ ) {
                index = p->index[j];
                m = ( header.mass[pt] < 1e-5 ) ? Particle[pt].m[index] : header.mass[pt];
                Particle1[pt].m[i] += m;
            }
            /*
            if ( 0 != Particle1[pt].m[i] && 1==pt ) {
                fprintf( stdout, "num=%i, dm=%f, m=%f\n", p->num, header.mass[pt], Particle1[pt].m[i] );
            }
            */
        }
    }
    fputs( sep_str, stdout );
}

void Generate_Header1() {
    long i, pt, num;;
    header1 = header;
    for ( pt=0; pt<2; pt++ ) {
        num = 0;
        for ( i=0; i<BINS3; i++ ) {
            if (Particle1[pt].m[i] > 1e-5)
                num++;
        }
        header1.npart[pt] = num;
        header1.npartTotal[pt] = num;
        header1.npartTotalHighWord[pt] = 0;
        header1.mass[pt] = 0;
    }
    for ( pt=2; pt<6; pt++ ) {
        header1.npart[pt] = 0;
        header1.npartTotal[pt] = 0;
        header1.mass[pt] = 0;
        header1.npartTotalHighWord[pt] = 0;
    }
    header1.num_files = 1;
}

void Allocate_Memory() {
    int pt;
    fputs( sep_str, stdout );
    fputs( "Change_NpartAllocate: Memory ...\n", stdout );
    for ( pt=0; pt<2; pt++ ) {
        Particle1[pt].pos = ( float* ) malloc( sizeof( float ) * header1.npart[pt] * 3 );
        Particle1[pt].vel = ( float* ) malloc( sizeof( float ) * header1.npart[pt] * 3 );
        Particle1[pt].id = ( long* ) malloc( sizeof( long ) * header1.npart[pt] );
        if ( 0 == pt )
            Particle1[pt].u = ( float* ) malloc( sizeof( float ) * header1.npart[pt] );
    }
    fputs( sep_str, stdout );
}

void Generate_Particle1() {
    long pt, i, j, index,id;
    float m;
    struct Particle_Line_Struct *p;
    fputs( sep_str, stdout );
    fputs( "Generate Particle1 Data ...\n", stdout );
    for ( pt=0; pt<2; pt++ ) {
        id = 0;
        Particle1[pt].num = header1.npart[pt];
        for ( i=0; i<BINS3; i++ ) {
            p = &( Particle_Line[pt][i] );
            if ( Particle1[pt].m[i] != 0 ) {
                Particle1[pt].id[id] = id;
                Particle1[pt].vel[id*3+0] = 0;
                Particle1[pt].vel[id*3+1] = 0;
                Particle1[pt].vel[id*3+2] = 0;
                Particle1[pt].pos[id*3+0] = 0;
                Particle1[pt].pos[id*3+1] = 0;
                Particle1[pt].pos[id*3+2] = 0;
                if ( 0 == pt ) {
                    Particle1[pt].u[id] = 0;
                }
                for ( j=0; j<p->num; j++ ) {
                    index = p->index[j];
                    m = ( header.mass[pt] < 1e-5 ) ? Particle[pt].m[index] : header.mass[pt];

                    Particle1[pt].pos[id*3+0] +=
                            Particle[pt].pos[index*3+0] * m;
                    Particle1[pt].pos[id*3+1] +=
                            Particle[pt].pos[index*3+1] * m;
                    Particle1[pt].pos[id*3+2] +=
                            Particle[pt].pos[index*3+2] * m;

                    Particle1[pt].vel[id*3+0] +=
                            Particle[pt].vel[index*3+0] * m;
                    Particle1[pt].vel[id*3+1] +=
                            Particle[pt].vel[index*3+1] * m;
                    Particle1[pt].vel[id*3+2] +=
                            Particle[pt].vel[index*3+2] * m;

                    if ( 0 == pt )
                        Particle1[pt].u[id] +=
                            Particle[pt].u[index] * m;
                }
                Particle1[pt].pos[id*3+0] /= Particle1[pt].m[i];
                Particle1[pt].pos[id*3+1] /= Particle1[pt].m[i];
                Particle1[pt].pos[id*3+2] /= Particle1[pt].m[i];

                Particle1[pt].vel[id*3+0] /= Particle1[pt].m[i];
                Particle1[pt].vel[id*3+1] /= Particle1[pt].m[i];
                Particle1[pt].vel[id*3+2] /= Particle1[pt].m[i];

                if ( 0 == pt )
                    Particle1[pt].u[id] /= Particle1[pt].m[i];
                id++;
                //fprintf( stdout, "id = %li\n", id );
            }
        }
    }
    fputs( sep_str, stdout );
}

void Update_Mass() {
    float *m, *tmp;
    long pt, i, id;
    fputs( sep_str, stdout );
    fputs( "Update Mass ...\n", stdout );
    FILE *fd;
    for ( pt=0; pt<2; pt++ ) {
        id = 0;
        m = ( float * ) malloc( sizeof( float ) * header1.npart[pt] );
        for ( i=0; i<BINS3; i++ ) {
            if ( Particle1[pt].m[i] != 0 ) {
                m[id] = Particle1[pt].m[i];
                id ++;
            }
        }
        tmp = Particle1[pt].m;
        Particle1[pt].m = m;
        free( tmp );
    }
    fputs( sep_str, stdout );
}

void Free_Memory() {
    long i, pt;
    fputs( sep_str, stdout );
    fputs( "Change_Npart: Free Memory ...\n", stdout );
    for ( pt=0; pt<2; pt++ ) {
        free( Particle1[pt].pos );
        free( Particle1[pt].vel );
        free( Particle1[pt].m );
        free( Particle1[pt].id );
        if ( 0 == pt )
            free( Particle1[pt].u );
        for ( i=0; i<BINS3; i++ ) {
            free( Particle_Line[pt][i].index );
        }
    }
    fputs( sep_str, stdout );
}

int Compare_For_Test_Data( const void *a, const void *b ) {
    return ( *(float*)b - *(float*)a > 0) ? 1 : 0;
}

void Test_Data() {
    long i, j, pt, num;
    float *m, tmp;
    fputs( sep_str, stdout );
    fprintf( stdout, "Test Data ...\n" );
    num = 10;
    for ( pt=0; pt<2; pt++ ) {
        m = ( float* ) malloc( sizeof( float ) * Particle1[pt].num );
        memcpy( m, Particle1[pt].m, sizeof( float ) * Particle1[pt].num );
        qsort( (void*)m, Particle1[pt].num, sizeof( float ), &Compare_For_Test_Data );
        fprintf( stdout, "Particle %i Mass:\n", pt );
        for ( i=0; i<num; i++ ) {
            fprintf( stdout, "%f\n", m[i] );
        }
        free( m );

    }
    fputs( sep_str, stdout );
}

void Change_Npart() {
    Construct_Particle_Line();
    Generate_Mass();
    Show_Header( header );
    Generate_Header1();
    Show_Header( header1 );
    Allocate_Memory();
    Generate_Particle1();
    Update_Mass();
    Write_File( Out_File, header1, Particle1 );
    Test_Data();
    Free_Memory();
}

