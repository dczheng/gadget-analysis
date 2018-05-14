void read_snapshot();
void free_memory();
void group_analysis();
void read_group();
void free_group();
void endrun( int ierr );
void read_parameters();
void set_units();

void slice();
void make_slice_img( int pt );
void write_img( char *fn );
void init_img();
void free_img();

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
void fof_save_groups();
int ngb_fof(  double *searchcenter, double h );

void conv( double ds );

