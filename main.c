#include "allvars.h"

void end_run( int ierr ) {
    fprintf( stderr, "EXIT CODE: %i\n", ierr );
    exit( ierr );
}

void read_para() {
    FILE *fd;
    int i;
    char *s;
    char line[ MAX_PARA_FILE_LINE_LEN ],
         name[ MAX_PARA_FILE_LINE_LEN ],
         data[ MAX_PARA_FILE_LINE_LEN ];
    fputs( sep_str, stdout );
    fprintf( stdout, "read Parameter... \n" );
    fd = fopen( para_file, "r" );
    if ( NULL == fd ) {
        fprintf( stderr, "Faile to Open Parameter file %s\n", para_file );
        end_run( 1 );
    }
    fgets( line, MAX_PARA_FILE_LINE_LEN, fd );
    slice_index_num = -1;
    while( !feof( fd ) ) {
        if ( ( line[0] != '#' ) && ( line[0] != '\n' ) ) {
            sscanf( line, "%s %s", name, data );
            if ( !strcmp( "FILE_PREFIX", name ) ) {
                strcpy( file_prefix, data );
                fprintf( stdout, "file_prefix: %s\n",  file_prefix );
            }
            if ( !strcmp( "GROUP_DIR", name ) ) {
                strcpy( group_dir, data );
                if ( group_dir[strlen(group_dir)-1] != '/' )
                    strcat( group_dir, "/" );
                fprintf( stdout, "group_dir: %s\n",  group_dir );
            }
            if ( !strcmp( "OUT_FILE", name ) ) {
                strcpy( out_file, data );
                fprintf( stdout, "out_file: %s\n",  out_file );
            }
            if ( !strcmp( "OUT_PICTURE_PREFIX", name ) ) {
                strcpy( out_picture_prefix, data );
                fprintf( stdout, "out_picture_prefix: %s\n",  out_picture_prefix );
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
            if ( !strcmp( "SCALAR_UNIT", name ) ) {
                scalar_unit = atof( data );
                fprintf( stdout, "scalar_unit: %e\n", scalar_unit );
            }
            if ( !strcmp( "3D_BOX", name ) ) {
                s = strtok( line, " " );
                fprintf( stdout, "box: " );
                for ( i=0; i<3; i++ ) {
                    if ( NULL == s ){
                        fprintf( stdout, "too few parameters for box\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    box[i] = atoi( s );
                    fprintf( stdout, "%i ", box[i] );
                }
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "3D_AL", name ) ) {
                s = strtok( line, " " );
                fprintf( stdout, "al: " );
                for ( i=0; i<3; i++ ) {
                    if ( NULL == s ){
                        fprintf( stdout, "too few parameters for al\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    al[i] = atof( s );
                    fprintf( stdout, "%.2f ", al[i] );
                }
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "3D_AZ", name ) ) {
                s = strtok( line, " " );
                fprintf( stdout, "az: " );
                for ( i=0; i<3; i++ ) {
                    if ( NULL == s ){
                        fprintf( stdout, "too few parameters for az\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    az[i] = atof( s );
                    fprintf( stdout, "%.2f ", az[i] );
                }
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "3D_CORNER1", name ) ) {
                s = strtok( line, " " );
                fprintf( stdout, "corner1: " );
                for ( i=0; i<3; i++ ) {
                    if ( NULL == s ){
                        fprintf( stdout, "too few parameters for corner1\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    corner1[i] = atof( s );
                    fprintf( stdout, "%.2f ", corner1[i] );
                }
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "3D_CORNER2", name ) ) {
                s = strtok( line, " " );
                fprintf( stdout, "corner2: " );
                for ( i=0; i<3; i++ ) {
                    if ( NULL == s ){
                        fprintf( stdout, "too few parameters for corner2\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    corner2[i] = atof( s );
                    fprintf( stdout, "%.2f ", corner2[i] );
                }
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "SLICE_INDEX_NUM", name ) ) {
                slice_index_num = atoi( data );
                fprintf( stdout, "slice_index_num: %i\n", slice_index_num );
                if ( slice_index_num !=0 ) {
                    slice_index = ( int* ) malloc( sizeof( int ) * slice_index_num );
                    memset( slice_index, -1, sizeof( int ) * slice_index_num );
                }
            }
            if ( !strcmp( "SLICE_INDEX", name ) ) {
                if ( -1 == slice_index_num ) {
                    fprintf( stderr, "SLICE_INDEX_NUM must appear before SLICE_INDEX!\n" );
                    end_run( 2 );
                }
                if ( slice_index_num !=0 ) {
                    s = strtok( line, " " );
                    fprintf( stdout, "slice_index: " );
                    i=0;
                    while ( s=strtok( NULL, " " ) ){
                        slice_index[i] = atoi( s );
                        fprintf( stdout, "%i ", slice_index[i] );
                        i++;
                    }
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
    strcpy( para_file, argv[1] );
    read_para();
    read_all_data();
    //group_analysis();
    //plot_scalar( 0, IO_MAG );
    //plot_scalar( 0, IO_MAG );
    magnetic_field_analysis();
    //density_analysis();
    //plot_position( 1 );
    //plot_3d_position( 4 );
    //plot_3d_multi( 3 );
    //plot_3d_scalar( 0, IO_ELEC );
    //velocity_analysis();
    free_all_memory();
    if ( slice_index_num >0 )
        free( slice_index );
}

