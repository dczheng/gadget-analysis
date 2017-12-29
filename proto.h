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
void init_projection();
void print_log( char *log );
void init_plot();
void free_plot();
void plot_slice();
void plot_syn_slice();
void analysis();
double kernel( double q );
void init_kernel_matrix();
void free_kernel_matrix();
long get_particle_offset( int pt );
long get_particle_num( int pt );

void tree_build( int pt );
void tree_free();


double H_a( double a );
double comoving_distance( double a );
double angular_distance( double a );
double luminosity_distance( double a );

void fof( int pt );
void fof_free();
int ngb_fof(  double *searchcenter, double h, int pt );
