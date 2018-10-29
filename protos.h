void read_snapshot();
void free_particle_memory();
void group_analysis();
void read_group();
void free_group();
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

void check_group_flag();

void test_radio();
void gas_state();
void gas_temperature();
void gas_density();
void mass_function();
void total_radio_spectrum();
void compute_temperature();

void init_tab_F();
void free_tab_F();
void test_F();
double particle_radio( double nu, long i );
double radio( double (*f)(double, void*), double *params,
        double B, double nu, double pmin, double pmax );

double qtrap( double (*func)( double, void* ),
        void *params, double a, double b, double err );


void powerspec();
double PowerSpec_Efstathiou(double k);

double sigma( sigma_struct ss );

double top_hat_filter( double k, void *params );
void top_hat_filter_k_limit( double *k0, double *k1, void *params );
void test_ps();
void test_sigma();

double PS_dndM( double a, double M );

void check_flag();

void compute_particle_radio();
void free_particle_radio();
