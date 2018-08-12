void read_snapshot();
void free_memory();
void group_analysis();
void read_group();
void free_group();
void endrun( int ierr );
void read_parameters( char *fn );
void set_units();
void init_sep_str();

void slice();
void make_slice_img( int pt );

void init_img();
void free_img();
void write_img( char *fn, char *s, int mode );
#define write_img1( fn, s )  write_img( fn, s, 0 )
#define write_img2( fn, s )  write_img( fn, s, 1 )


void create_dir( char *s );
double second();

void analysis();
double kernel( double q );
void init_kernel_matrix();
void free_kernel_matrix();

long get_particle_offset( int pt );
long get_particle_num( int pt );
double get_B( long i );

void tree_build();
void tree_free();


double hubble_function( double a );
double comoving_distance( double a );
double angular_distance( double a );
double luminosity_distance( double a );
void compute_cosmo_quantities();

void fof();
void fof_free();
void fof_save();
void fof_read();
int ngb_fof(  double *searchcenter, double h );

void conv( double ds );
void init_conv_kernel();
void free_conv_kernel();

void check_flags();
void task_sync_test();
