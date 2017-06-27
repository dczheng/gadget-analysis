#include "allvars.h"

void end_run( int ierr ) {
    fprintf( stderr, "EXIT CODE: %i\n", ierr );
    exit( ierr );
}

void read_para() {
    FILE *fd;
    char line[ MAX_PARA_FILE_LINE_LEN ],
         name[ MAX_PARA_FILE_LINE_LEN ],
         data[ MAX_PARA_FILE_LINE_LEN ];
    fputs( sep_str, stdout );
    fprintf( stdout, "read Parameter... \n" );
    fd = fopen( Para_file, "r" );
    if ( NULL == fd ) {
        fprintf( stderr, "Faile to Open Parameter file %s\n", Para_file );
        end_run( 1 );
    }
    fgets( line, MAX_PARA_FILE_LINE_LEN, fd );
    slice_index_num = -1;
    while( !feof( fd ) ) {
        if ( ( line[0] != '#' ) && ( line[0] != '\n' ) ) {
            sscanf( line, "%s %s", name, data );
            if ( !strcmp( "FILE_PREFIX", name ) ) {
                strcpy( file_Prefix, data );
                fprintf( stdout, "file_Prefix: %s\n",  file_Prefix );
            }
            if ( !strcmp( "OUT_FILE", name ) ) {
                strcpy( Out_file, data );
                fprintf( stdout, "Out_file: %s\n",  Out_file );
            }
            if ( !strcmp( "OUT_PICTURE_PREFIX", name ) ) {
                strcpy( Out_Picture_Prefix, data );
                fprintf( stdout, "Out_Picture_Prefix: %s\n",  Out_Picture_Prefix );
            }
            if ( !strcmp( "NUM_FILES", name ) ) {
                Num_files = atoi( data );
                fprintf( stdout, "Num_files: %i\n",  Num_files );
            }
            if ( !strcmp( "PIC_XSIZE", name ) ) {
                pic_xsize = atoi( data );
                fprintf( stdout, "pic_xsize: %i\n",  pic_xsize );
            }
            if ( !strcmp( "PIC_YSIZE", name ) ) {
                pic_ysize = atoi( data );
                fprintf( stdout, "pic_ysize: %i\n",  pic_ysize );
            }
            if ( !strcmp( "SLICE_NUM", name ) ) {
                slice_num = atoi( data );
                fprintf( stdout, "slice_num: %i\n", slice_num );
            }
            if ( !strcmp( "REDSHIFT", name ) ) {
                redshift = atof( data );
                fprintf( stdout, "reshift: %.2f\n", redshift );
            }
            if ( !strcmp( "SLICE_INDEX_NUM", name ) ) {
                slice_index_num = atoi( data );
                fprintf( stdout, "slice_index_num: %i\n", slice_index_num );
                slice_index = ( int* ) malloc( sizeof( int ) * slice_index_num );
                memset( slice_index, -1, sizeof( int ) * slice_index_num );
            }
            if ( !strcmp( "SLICE_INDEX", name ) ) {
                if ( -1 == slice_index_num ) {
                    fprintf( stderr, "SLICE_INDEX_NUM must appear before SLICE_INDEX!\n" );
                    end_run( 2 );
                }
                int i;
                char *s;
                s = strtok( line, " " );
                fprintf( stdout, "slice_index: " );
                i=0;
                while ( s=strtok( NULL, " " ) ){
                    slice_index[i] = atoi( s );
                    fprintf( stdout, "%i ", slice_index[i] );
                    i++;
                }
                fprintf( stdout, "\n" );
            }
        }
        fgets( line, MAX_PARA_FILE_LINE_LEN, fd );
    }
    fputs( sep_str, stdout );
    fclose( fd );
}

void init_sep_str() {
    memset( sep_str, '-', SEP_LEN-1 );
    sep_str[ SEP_LEN-1 ] = '\n';
}

int main( int argc, char *argv[] ){
    int i;
    if ( argc < 2 ) {
        fprintf( stderr, "Parameter file is required on command line!\n " );
        end_run( 1 );
    }
    init_sep_str();
    strcpy( Para_file, argv[1] );
    read_para();
    read_all_data();
    plot_scalar( 0, IO_RHO );
    //plot_scalar( 0, IO_MAG );
    //Analysis_Magnetic_Field();
    //Plot_2D_Point( 0 );
    free_all_memory();
    free( slice_index );
}

