void read_snapshot();
void free_particle_memory();
void group_analysis();
void read_group();
void free_group();
void read_parameters( char *fn );
void set_units();

void init_img();
void free_img();
void reset_img();
void write_img( char *fn, char *s, int mode );
#define write_img1( fn, s )  write_img( fn, s, 0 )
#define write_img2( fn, s )  write_img( fn, s, 1 )

void create_dir( char *s );
double second();

void analysis();
double kernel( double q );
void init_kernel_matrix();
void free_kernel_matrix();

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

void check_group_flag();

void test_radio();
void phase();

void slice_init();
void make_slice_img( int pt, double *data );

void temperature_slice();
void radio_slice();
void density_slice();
void mach_slice();
void mag_slice();
void cren_slice();
void mass_function();
void total_radio_spectrum();
void compute_temperature();

void init_tab_F();
void free_tab_F();
void test_F();
void init_compute_F();
double particle_radio( double nu, long i );

double radio_inte( double p, void *params );
double radio( double (*f)(double, void*), double *params,
        double B, double nu, double pmin, double pmax, double err );

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
void B_Pdf();
void pre_proc();

void do_sync( char *s );
void do_sync_local( char *s );
void do_sync_master( char *s );
void task_sync_test( char *s );

void test_cos();

void cre_pressure_pdf();

void part_info();

void corr_Tdiff_dens();
void corr_dens();

void field_to_grid( double *data, double *grid, double *num, long N, int flag );
