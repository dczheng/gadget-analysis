#include "allvars.h"

char buf[100], *bname;

void init_sep_str() {
    memset( sep_str, '-', SEP_LEN-2 );
    sep_str[ SEP_LEN-2 ] = '\n';
    sep_str[ SEP_LEN-1 ] = '\0';

    memset( sep_str0, '-', SEP_LEN0-2 );
    sep_str0[ SEP_LEN0-2 ] = '\n';
    sep_str0[ SEP_LEN0-1 ] = '\0';
}

void global_init() {
    init_sep_str();
    ms.max_mem = ms.mem = ms.nn = 0;
    srand( time(NULL) );
    inte_ws = gsl_integration_workspace_alloc( GSL_INTE_WS_LEN );
}

void global_free() {
    gsl_integration_workspace_free( inte_ws );
}

void mpi_comms_test() {

    int a, *arr, arrN, i, windisp, k;
    MPI_Aint winsize;
    MPI_Win MpiWin;

    do_sync( "" );

    put_sep;
    writelog( "mpi comms test ...\n" );
    arrN = 10;

    if ( ThisTask_Local == 0 ) {
        a = ThisTask * 20;
    }

    MPI_Bcast( &a, 1, MPI_INT, 0, MpiComm_Local );
    do_sync( "" );

    for ( k=0; k<NTask/NTask_Local; k++ ) {
        if ( ( k * NTask_Local <= ThisTask ) && ( ThisTask < (k+1) * NTask_Local ) )
            printf( "ThisTask: %i, a: %i\n", ThisTask, a );
        do_sync( "" );
    }

    if ( ThisTask_Master == 0 )
        a = 99999;

    if ( ThisTask_Local == 0 )
        MPI_Bcast( &a, 1, MPI_INT, 0, MpiComm_Master );

    do_sync( "" );

    if ( ThisTask_Local == 0 )
        printf( "ThisTask: %i, ThisTask_Master: %i, a: %i\n",
                ThisTask, ThisTask_Master, a );

    do_sync( "" );

    MPI_Win_allocate_shared( arrN * sizeof(int), sizeof(int),
            MPI_INFO_NULL, MpiComm_Local, &arr, &MpiWin );

    MPI_Win_fence( 0, MpiWin );

    if ( ThisTask_Local == 0 ) {
        for( i=0; i<arrN; i++ ) {
            arr[i] = ThisTask_Master * arrN + i;
        }
    }
    else {
        MPI_Win_shared_query( MpiWin, 0, &winsize, &windisp, &arr );
    }

    MPI_Win_fence( 0, MpiWin );

    for( k=0; k<NTask; k++ ) {
        for( i=0; i<arrN; i++ ) {
            if ( i % NTask_Local == ThisTask_Local ) {
                if ( ThisTask == k )
                printf( "ThisTask: %i, ThisTask_Local: %i, get: %i\n",
                        ThisTask, ThisTask_Local, arr[i] );
            }
        }
        do_sync( "" );

    }

    do_sync( "" );

    endrun( 20181207 );
}

void create_mpi_comms() {

    int range[3];
    MPI_Group g, g0;
    MPI_Comm  TmpComm;


    writelog( "create mpi comms ...\n" );

    NTask_Local = NumThreadsPerSnapshot;
    NTask_Master = NTask / NTask_Local;

    MasterTask = ThisTask / NTask_Local * NTask_Local;

    range[0] = 0;
    range[1] = NTask-1;
    range[2] = NTask_Local;

    MPI_Comm_group( MPI_COMM_WORLD, &g0 );
    MPI_Group_range_incl( g0, 1, &range, &g );

    MPI_Comm_create( MPI_COMM_WORLD, g, &MpiComm_Master );

    MPI_Group_free( &g0 );
    MPI_Group_free( &g );


    ThisTask_Master = -99999;

    if ( ThisTask == MasterTask ) {
        MPI_Comm_rank( MpiComm_Master, &ThisTask_Master );
        /*
        printf( "ThisTask: %i, ThisTask_Master: %i\n",
            ThisTask, ThisTask_Master );
            */
    }

    MPI_Comm_split( MPI_COMM_WORLD, MasterTask, ThisTask, &TmpComm );
    MPI_Comm_split_type( TmpComm, MPI_COMM_TYPE_SHARED, ThisTask,
            MPI_INFO_NULL, &MpiComm_Local );
    MPI_Comm_rank( MpiComm_Local, &ThisTask_Local );
    MPI_Comm_free( &TmpComm );

    /*
    printf( "ThisTask: %i, ThisTask_Local: %i, ThisTask_SharedMem: %i\n",
            ThisTask, ThisTask_Local, ThisTask_SharedMem );
            */

    //mpi_comms_test();

    do_sync( "create comm" );

}

void free_comms() {

    writelog( "free mpi comms ...\n" );
    MPI_Comm_free( &MpiComm_Local );
    if ( ThisTask_Local == 0 )
        MPI_Comm_free( &MpiComm_Master );

}

void merge_log_file() {

    int i, j, status;
    FILE *fd,*fd1;
    char buf2[200], t[2][10];

    t[0][0]  = '\0';
    sprintf( t[1], "-usedmem");

    printf( "merge log file ...\n" );
    for ( i=0; i<2; i++ ) {
        sprintf( buf, "./gadget-analysis.log/%s%s.log", bname, t[i] );
        //printf( "%s\n", buf );
        fd = fopen( buf, "w" );
        if ( NULL == fd ) {
            printf( "Failed to open file %s\n", buf );
            endrun( 20181213 );
        }

        for( j=0; j<NTask; j++ ) {
            sprintf( buf, "./gadget-analysis.log/%s%s-%03i.log", bname, t[i], j );
       // printf( "%s\n", buf );

            fd1 = fopen( buf, "r" );
            if ( NULL == fd1 ) {
                printf( "Failed to open file %s\n", buf );
                endrun( 20181213 );
            }

            fprintf( fd, "%s", sep_str );
            fprintf( fd, "Task: %i\n", j );
            fprintf( fd, "%s", sep_str );

            while( !feof( fd1 ) ) {
                fgets( buf2, 200, fd1 );
                fprintf( fd, buf2 );
            }

            fclose( fd1 );

            status = remove( buf );

            if ( status ) {
                printf( "Failed to remove file %s\n", buf );
                endrun( 20181213 );
            }

        }

        fclose( fd );
    }

    //endrun( 20181213 );

}

int main( int argc, char *argv[] ){

    time_t time1, time2;
    long dtime;
    struct tm *tb;
    int provided;

    //MPI_Init( &argc, &argv );
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &provided );
    //MPI_Init_thread( &argc, &argv, MPI_THREAD_SERIALIZED, &provided );
    MPI_Comm_rank( MPI_COMM_WORLD, &ThisTask );
    MPI_Comm_size( MPI_COMM_WORLD, &NTask );

    if ( ThisTask == 0 ) {

        if ( argc < 2 )
            endruns( "Parameters file is required on command line! " );

        if ( access( argv[1], 0 ) == -1 )
            endrun0( "Parameters file `%s` is invalid!\n", argv[1] );

    }

    if ( argc >= 3 )
        NumThreadsPerSnapshot = atoi( argv[2] );
    else
        NumThreadsPerSnapshot = 1;

    if ( NTask % NumThreadsPerSnapshot != 0 ) {
        printf( "NTask `%i` must be NumThreadsPerSnapshot `%i` * n\n",
                NTask, NumThreadsPerSnapshot );
        endrun( 20181205 );
    }

    time1 = time( NULL );
    tb = localtime( &time1 );

    global_init();

    bname = basename( argv[1] );
    //printf( "argv[1]: %s, bname: %s\n", argv[1], bname );

    /************************init log**************************************/
    if ( ThisTask == 0 ){
        printf( "%s", sep_str );
        if ( access( "./gadget-analysis.log/", 0 ) == -1 ) {
            printf( "create directory `./gadget-analysis.log/` by task 0\n" );

            if ( mkdir( "./gadget-analysis.log/", 0755 ) == -1 )
                endruns( "failed create directory `./gadget-analysis.log/`" );
        }
    }

    MPI_Barrier( MPI_COMM_WORLD );

    /*
    if ( ThisTask == 0 )
    merge_log_file();
    */

    sprintf( buf, "./gadget-analysis.log/%s-%03d.log", bname, ThisTask );
    LogFileFd = fopen( buf, "w" );

    if ( NULL == LogFileFd )
        endrun0( "Failed open file %s\n", buf );

    sprintf( buf, "./gadget-analysis.log/%s-usedmem-%03d.log", bname, ThisTask );
    UsedMemFileFd = fopen( buf, "w" );

    if ( NULL == UsedMemFileFd )
        endrun0( "Failed open file %s\n", buf );

    sprintf( OutputDir, "./output_%s/", bname );

    create_mpi_comms();
    put_sep;

    writelog( "Start At: %s", asctime(tb) );

#ifdef ZDEBUG
    writelog( "Assign `SIGSEGV` to signal_hander function.\n" );
    signal( SIGSEGV, signal_hander );
    init_sig();
#endif

    ToolsPath = getenv( "GADGET_TOOLS" );
    if ( strcmp( ToolsPath, "" ) == 0 ){
        writelog( "Please set `GADGET_TOOLS` evironment variable.\n" );
        put_sep;
        endrun(20181107);
    }
    writelog( "GADGET_TOOLS: %s\n", ToolsPath );
    writelog( "NumThreadsPerSnapshot: %i\n", NumThreadsPerSnapshot );
    writelog( "Parameters file: %s\n", argv[1] );


    /******************read***********************/

    read_parameters( argv[1] );
    check_flag();

    if ( All.ParallelIO <= 0 ) {
        All.ParallelIO = NTask_Master;
    }

    IOGroups = NTask_Master / All.ParallelIO;
    if ( NTask_Master % All.ParallelIO )
        IOGroups ++;

    //test_sigma();
    //
    //test_radio();

    //test_cos();

    read_snapshot();
    do_sync( "read data" );

    /*
    int i;

    for ( i=0; i<N_Gas; i++ ){
//        SphP[i].CRE_Alpha = 2.1;
//        SphP[i].CRE_C *= 1e10;

        SphP[i].CRE_qmax = 1e8;
        //SphP[i].CRE_qmin = 1;
        //SphP[i].CRE_qmax *= 3.1;
    }
    */

    pre_proc();
    do_sync( "pre process data" );

    /******************read***********************/

    set_units();

    compute_cosmo_quantities();

    //test_ps();

    /******************analysis***********************/

    analysis();
    /******************analysis***********************/
    do_sync( "analysis" );

    free_particle_memory();
    do_sync( "" );

    free_comms();

    time2 = time( NULL );
    tb = localtime( &time2 );

    writelog( "End At: %s", asctime(tb) );

    dtime = (long) ( difftime( time2, time1 ) );
    writelog( "Total Time %li sec.\n", dtime );
    put_sep;

    fclose( LogFileFd );
    fclose( UsedMemFileFd );

    global_free();

    if ( ThisTask == 0 )
        merge_log_file();

    MPI_Finalize();
    return 0;
}

