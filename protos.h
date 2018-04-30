void read_snapshot();
void free_memory();
void group_analysis();
void read_group();
void free_group();
void signal_hander( int sig );
void endrun( int ierr );
void read_parameters();
void set_units();

void slice();
void make_slice_img( int pt );
void write_img( char *fn, int mode);
void create_dir( char *s );

void analysis();
double kernel( double q );
void init_kernel_matrix();
void free_kernel_matrix();

long get_particle_offset( int pt );
long get_particle_num( int pt );
double get_B( long i );

void tree_build( int pt );
void tree_free();


double hubble_function( double a );
double comoving_distance( double a );
double angular_distance( double a );
double luminosity_distance( double a );
void compute_cosmo_quantities();

void fof( int pt );
void fof_free();
void fof_save_groups();
int ngb_fof(  double *searchcenter, double h, int pt );
