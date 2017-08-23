#include "allvars.h"

void end_run( int ierr ) {
if ( this_task == 0 )
    fprintf( stderr, "EXIT CODE: %i\n", ierr );
    MPI_Abort( MPI_COMM_WORLD, ierr );
    exit( ierr );
}

void read_para() {
    FILE *fd;
    int i;
    char *s;
    char line[ MAX_PARA_FILE_LINE_LEN ],
         name[ MAX_PARA_FILE_LINE_LEN ],
         data[ MAX_PARA_FILE_LINE_LEN ];
if ( this_task == 0 )
    fputs( sep_str, stdout );
if ( this_task == 0 )
    fprintf( stdout, "read Parameter... \n" );
    fd = fopen( para_file, "r" );
    if ( NULL == fd ) {
if ( this_task == 0 )
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
if ( this_task == 0 )
                fprintf( stdout, "file_prefix: %s\n",  file_prefix );
            }
            if ( !strcmp( "GROUP_DIR", name ) ) {
                strcpy( group_dir, data );
                if ( group_dir[strlen(group_dir)-1] != '/' )
                    strcat( group_dir, "/" );
if ( this_task == 0 )
                fprintf( stdout, "group_dir: %s\n",  group_dir );
            }
            if ( !strcmp( "OUT_FILE", name ) ) {
                strcpy( out_file, data );
if ( this_task == 0 )
                fprintf( stdout, "out_file: %s\n",  out_file );
            }
            if ( !strcmp( "OUT_PICTURE_PREFIX", name ) ) {
                strcpy( out_picture_prefix, data );
if ( this_task == 0 )
                fprintf( stdout, "out_picture_prefix: %s\n",  out_picture_prefix );
            }
            if ( !strcmp( "FLAG_ARROW", name ) ) {
                flag_arrow = atoi( data );
if ( this_task == 0 )
                fprintf( stdout, "flag_arrow: %i\n",  flag_arrow );
            }
            if ( !strcmp( "ARROW_X", name ) ) {
                arrow_x = atoi( data );
if ( this_task == 0 )
                fprintf( stdout, "arrow_x: %i\n",  arrow_x );
            }
            if ( !strcmp( "ARROW_Y", name ) ) {
                arrow_y = atoi( data );
if ( this_task == 0 )
                fprintf( stdout, "arrow_y: %i\n",  arrow_y );
            }
            if ( !strcmp( "NUM_FILES", name ) ) {
                Num_files = atoi( data );
if ( this_task == 0 )
                fprintf( stdout, "Num_files: %i\n",  Num_files );
            }
            if ( !strcmp( "PROJECTION_MODE", name ) ) {
                proj_mode = atoi( data );
if ( this_task == 0 )
                fprintf( stdout, "projection mode: %i\n",  proj_mode );
            }
            if ( !strcmp( "OUT_PIC_DATA", name ) ) {
                out_pic_data = atoi( data );
if ( this_task == 0 )
                fprintf( stdout, "out_pic_data: %i\n",  out_pic_data );
            }
            if ( !strcmp( "PIC_XSIZE", name ) ) {
                pic_xsize = atoi( data );
if ( this_task == 0 )
                fprintf( stdout, "pic_xsize: %i\n",  pic_xsize );
            }
            if ( !strcmp( "PIC_YSIZE", name ) ) {
                pic_ysize = atoi( data );
if ( this_task == 0 )
                fprintf( stdout, "pic_ysize: %i\n",  pic_ysize );
            }
            if ( !strcmp( "SLICE_NUM", name ) ) {
                slice_num = atoi( data );
if ( this_task == 0 )
                fprintf( stdout, "slice_num: %i\n", slice_num );
            }
            if ( !strcmp( "REDSHIFT", name ) ) {
                redshift = atof( data );
if ( this_task == 0 )
                fprintf( stdout, "reshift: %.2f\n", redshift );
            }
            if ( !strcmp( "SCALAR_UNIT", name ) ) {
                scalar_unit = atof( data );
if ( this_task == 0 )
                fprintf( stdout, "scalar_unit: %e\n", scalar_unit );
            }
            if ( !strcmp( "3D_BOX", name ) ) {
                s = strtok( line, " " );
if ( this_task == 0 )
                fprintf( stdout, "box: " );
                for ( i=0; i<3; i++ ) {
                    if ( NULL == s ){
if ( this_task == 0 )
                        fprintf( stdout, "too few parameters for box\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    box[i] = atoi( s );
if ( this_task == 0 )
                    fprintf( stdout, "%i ", box[i] );
                }
if ( this_task == 0 )
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "3D_AL", name ) ) {
                s = strtok( line, " " );
if ( this_task == 0 )
                fprintf( stdout, "al: " );
                for ( i=0; i<3; i++ ) {
                    if ( NULL == s ){
if ( this_task == 0 )
                        fprintf( stdout, "too few parameters for al\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    al[i] = atof( s );
if ( this_task == 0 )
                    fprintf( stdout, "%.2f ", al[i] );
                }
if ( this_task == 0 )
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "3D_AZ", name ) ) {
                s = strtok( line, " " );
if ( this_task == 0 )
                fprintf( stdout, "az: " );
                for ( i=0; i<3; i++ ) {
                    if ( NULL == s ){
if ( this_task == 0 )
                        fprintf( stdout, "too few parameters for az\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    az[i] = atof( s );
if ( this_task == 0 )
                    fprintf( stdout, "%.2f ", az[i] );
                }
if ( this_task == 0 )
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "SLICE_CORNER1", name ) ) {
                s = strtok( line, " " );
if ( this_task == 0 )
                fprintf( stdout, "slice_corner1: " );
                for ( i=0; i<2; i++ ) {
                    if ( NULL == s ){
if ( this_task == 0 )
                        fprintf( stdout, "too few parameters for slice_corner1\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    slice_corner1[i] = atof( s );
if ( this_task == 0 )
                    fprintf( stdout, "%.2f ", slice_corner1[i] );
                }
if ( this_task == 0 )
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "SLICE_CORNER2", name ) ) {
                s = strtok( line, " " );
if ( this_task == 0 )
                fprintf( stdout, "sclie_corner2: " );
                for ( i=0; i<2; i++ ) {
                    if ( NULL == s ){
if ( this_task == 0 )
                        fprintf( stdout, "too few parameters for slice_corner2\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    slice_corner2[i] = atof( s );
if ( this_task == 0 )
                    fprintf( stdout, "%.2f ", slice_corner2[i] );
                }
if ( this_task == 0 )
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "3D_CORNER1", name ) ) {
                s = strtok( line, " " );
if ( this_task == 0 )
                fprintf( stdout, "corner1: " );
                for ( i=0; i<3; i++ ) {
                    if ( NULL == s ){
if ( this_task == 0 )
                        fprintf( stdout, "too few parameters for corner1\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    corner1[i] = atof( s );
if ( this_task == 0 )
                    fprintf( stdout, "%.2f ", corner1[i] );
                }
if ( this_task == 0 )
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "3D_CORNER2", name ) ) {
                s = strtok( line, " " );
if ( this_task == 0 )
                fprintf( stdout, "corner2: " );
                for ( i=0; i<3; i++ ) {
                    if ( NULL == s ){
if ( this_task == 0 )
                        fprintf( stdout, "too few parameters for corner2\n" );
                        end_run( 2 );
                    }
                    s = strtok( NULL, " " );
                    corner2[i] = atof( s );
if ( this_task == 0 )
                    fprintf( stdout, "%.2f ", corner2[i] );
                }
if ( this_task == 0 )
                fprintf( stdout, "\n" );
            }
            if ( !strcmp( "SLICE_INDEX_NUM", name ) ) {
                slice_index_num = atoi( data );
if ( this_task == 0 )
                fprintf( stdout, "slice_index_num: %i\n", slice_index_num );
                if ( slice_index_num !=0 ) {
                    slice_index = ( int* ) malloc( sizeof( int ) * slice_index_num );
                    memset( slice_index, -1, sizeof( int ) * slice_index_num );
                }
            }
            if ( !strcmp( "SLICE_INDEX", name ) ) {
                if ( -1 == slice_index_num ) {
if ( this_task == 0 )
                    fprintf( stderr, "SLICE_INDEX_NUM must appear before SLICE_INDEX!\n" );
                    end_run( 2 );
                }
                if ( slice_index_num !=0 ) {
                    s = strtok( line, " " );
if ( this_task == 0 )
                    fprintf( stdout, "slice_index: " );
                    i=0;
                    while ( s=strtok( NULL, " " ) ){
                        slice_index[i] = atoi( s );
if ( this_task == 0 )
                        fprintf( stdout, "%i ", slice_index[i] );
                        i++;
                    }
                }
if ( this_task == 0 )
                fprintf( stdout, "\n" );
            }
        }
        fgets( line, MAX_PARA_FILE_LINE_LEN, fd );
    }
if ( this_task == 0 )
    fputs( sep_str, stdout );
    fclose( fd );
}

void init_sep_str() {
    memset( sep_str, '-', SEP_LEN-1 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';
}

int main( int argc, char *argv[] ){
    int i;
    time_t time1, time2;
    struct tm *tb;
    char tmp[100];
    if ( argc < 2 ) {
if ( this_task == 0 )
        fprintf( stderr, "Parameter file is required on command line!\n " );
        end_run( 1 );
    }
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &this_task );
    MPI_Comm_size( MPI_COMM_WORLD, &task_num );
#ifdef DEBUG
    signal( SIGSEGV, signal_hander );
#endif
    time1 = time( NULL );
    tb = localtime( &time1 );
    if ( this_task == 0 )
        fprintf( stdout, "Start At: %s", asctime(tb) );
    init_sep_str();
    strcpy( para_file, argv[1] );
    read_para();
    read_all_data();
    //group_analysis();
    sprintf( tmp, "%s", out_picture_prefix );
    /*
    sprintf( out_picture_prefix, "%s/%s_%.2f", tmp, "mag", redshift );
    plot_slice( 0, IO_MAG );
    */
    /*
    sprintf( out_picture_prefix, "%s/%s_%.2f", tmp, "gas", redshift );
    plot_slice( 0, IO_MASS );
    */
    /*
    sprintf( out_picture_prefix, "%s/%s_%.2f", tmp, "dm", redshift );
    plot_slice( 1, IO_MASS );
    */
    /*
    if ( header.npart[4] != 0 ) {
        sprintf( out_picture_prefix, "%s/%s_%.2f", tmp, "star", redshift );
        plot_slice( 4, IO_MASS );
    }
    */
    /*
    sprintf( out_picture_prefix, "%s/%s_%.2f", tmp, "mn", redshift );
    plot_slice( 0, IO_MN );
    */
    sprintf( out_picture_prefix, "%s", tmp );
   //analysis_radio();
    //plot_slice( 0, IO_MAG );
    //magnetic_field_analysis();
    //density_analysis();
    //plot_position( 1 );
    plot_3d_position( 1 );
    //plot_3d_multi( 1 );
    //plot_3d_scalar( 0, IO_ELEC );
    //velocity_analysis();
    free_all_memory();
    if ( slice_index_num >0 )
        free( slice_index );
    MPI_Finalize();
    time2 = time( NULL );
    tb = localtime( &time2 );
    if ( this_task == 0 ){
        fprintf( stdout, "End At: %s", asctime(tb) );
        fprintf( stdout, "Total Time %i\n", time2-time1 );
    }
}

