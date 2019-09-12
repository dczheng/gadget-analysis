//%------>>>>>>file : conv.c
//%
double kernel_func( double x );
void init_conv_kernel( double ds );
void free_conv_kernel();
void conv( double ds );


//%------>>>>>>file : radio.c
//%
double F_integrand( double x, void *params );
void F0( double x, double *r, double *err );
void init_compute_F();
void F( double x, double *r, double *err );
double tab_F( double x );
void test_tab_F();
void output_F_x();
void output_tab_F();
void init_tab_F();
void free_tab_F();
double radio_inte( double p, void *params );
double radio( double (*f)( double, void* ), double *params,double B, double nu, double pmin, double pmax, double err );


//%------>>>>>>file : set_units.c
//%
void set_units();


//%------>>>>>>file : pre_proc.c
//%
long get_particle_num( int pt );
long get_particle_offset( int pt );
void find_id();
void construct_id_to_index();
int compare_pos( const void *a, const void *b );
int compare_sphp_pos( const void *a, const void *b );
void sort_particle_by_pos();
void test_sort();
void merge_particle( int pt );
void pre_proc();
void check_data( int err );


//%------>>>>>>file : part_radio.c
//%
double particle_df( double p, void *params );
double particle_radio2( double nu,  SphParticleData *part );
double particle_radio( double nu, long i );
void save_particle_radio();
int read_particle_radio();
void compute_particle_radio();
void free_particle_radio();
void d2PdVdv_qmax();
void output_radio_inte();
void test_part_radio();


//%------>>>>>>file : mf.c
//%
void mass_function();


//%------>>>>>>file : cre.c
//%
double beta_inte( double x, void *params );
double beta( double a, double b, double x );
void compute_cre_pressure();
void cre_pressure_pdf();


//%------>>>>>>file : tree.c
//%
void tree_allocate();
void tree_free();
void tree_build_single();
void tree_walk_recursive( long n, long sib, long father );
void tree_walk_test();
void tree_build();


//%------>>>>>>file : group_analysis.c
//%
int group_present( long index );
inline double get_group_size( struct group_properties *g );
void group_analysis();


//%------>>>>>>file : mymath.c
//%
double trapzd( double (*func)( double, void* ),void *params, double a, double b, int n );
double qtrap( double (*func)( double, void* ),void *params, double a, double b, double err );


//%------>>>>>>file : ngb.c
//%
int ngb_fof( double *searchcenter, double h );
int ngb( double *searchcenter, double h );


//%------>>>>>>file : part_info.c
//%
void part_info();


//%------>>>>>>file : total_radio_spectrum.c
//%
void total_radio_spectrum();


//%------>>>>>>file : smooth.c
//%
int smooth_flag();
void smooth();


//%------>>>>>>file : check_flag.c
//%
void check_flag();


//%------>>>>>>file : img.c
//%
void init_img();
void free_img();
void reset_img();
void write_img( char *fn, char *nstr );


//%------>>>>>>file : fof.c
//%
void fof_allocate( long N );
void fof_free();
int fof_compare_len( const void *a, const void *b );
int fof_compare_mass( const void *a, const void *b );
void fof_find_groups();
void fof_compute_group_properties();
void fof_test();
void fof_save();
void fof_read();
void fof();


//%------>>>>>>file : group_radio.c
//%
double group_luminosity( int nu_index, long index );
void group_flux( int nu_index, long index, double *flux, double *flux_nosr );
void group_spectrum();
double particle_f( SphParticleData *part, double p );
void group_electron_spectrum();
void group_spectrum_index();


//%------>>>>>>file : phase.c
//%
void phase();


//%------>>>>>>file : analysis.c
//%
void init_analysis();
void free_analysis();
void analysis();


//%------>>>>>>file : powerspec.c
//%
double powerspec_interp( double k );
void test_powerspec_interp();
void compute_sigma8();
void powerspec();


//%------>>>>>>file : grid.c
//%
void field_to_grid( double *data, double *grid, int pt, int Nmin, int flag );
void data_to_grid2d( double *data, double *grid, long Ndata, int NGrid, double L );


//%------>>>>>>file : group_temp.c
//%
void group_temp_profile();
void group_temp_stack();


//%------>>>>>>file : debug.c
//%
void signal_hander( int s );
void empty_sig_buf();
void init_sig();


//%------>>>>>>file : system.c
//%
void create_dir0( char *s );
void do_sync0( char *s, MPI_Comm comm );
void do_sync( char *s );
void do_sync_local( char *s );
void do_sync_master( char *s );
double second();
void task_sync_test( char *s );


//%------>>>>>>file : allvars.c
//%


//%------>>>>>>file : field.c
//%
void field_cren_T_dens();


//%------>>>>>>file : main.c
//%
void init_sep_str();
void global_init();
void global_free();
void mpi_comms_test();
void create_mpi_comms();
void free_comms();
void merge_log_file();
int main( int argc, char *argv[] );


//%------>>>>>>file : read_snapshot.c
//%
int blockpresent0( enum iofields blk, int pt );
int blockpresent( enum iofields blk, int pt );
int get_block_nbytes( enum iofields blk );
void get_block_dims( int pt, enum iofields blk, hsize_t (*dims)[2] );
void get_dataset_name( enum iofields blk, char *buf );
void get_hdf5_native_type( enum iofields blk, hid_t *hdf5_type );
void empty_buffer( enum iofields blk, int offset, int pt );
void read_header( char *fn );
void show_header( io_header header );
void write_header( char *fn, io_header header );
void free_particle_memory();
void read_snapshot_test();
void read_snapshot();


//%------>>>>>>file : temp.c
//%
void compute_temperature();


//%------>>>>>>file : group_plot.c
//%
int group_filed_present( enum group_fields blk );
void get_group_filed_name( enum group_fields blk, char *buf );
void group_plot();


//%------>>>>>>file : pdf.c
//%
void B_Pdf();
void dens_pdf();
void T_pdf();
void pdf2d_or_field2d( double *x, double *y, double *w, long num, char *dn,int flag, double *mm, int Nmin );
void cren_T_pdf();
void hsml_T_pdf();
void u_T_pdf();
void hsml_dens_pdf();


//%------>>>>>>file : gas_ratio.c
//%
void gas_ratio();


//%------>>>>>>file : kernel.c
//%
double kernel( double q );
void init_kernel_matrix();
void free_kernel_matrix();
inline void kernel_hinv(double h, double *hinv, double *hinv3, double *hinv4);
inline void kernel_main(double u, double hinv3, double hinv4,double *wk, double *dwk, int mode);


//%------>>>>>>file : correlation.c
//%
void set_global_vars();
void corr( double *x, double *y, double *c );
void get_corr1d( double *corr3d, double *corr1d );
void output_grid_slice( double *grid, int N, char *fn_prefix );
void put_dm_on_grid( double *grid );
void put_gas_on_grid( double *grid );
void put_dens_on_grid( double *grid, int Nmin );
void put_temp_on_grid( double *grid, int Nmin );
void corr_dm();
void corr_gas();
void corr_Tdiff_dens();
void pdf_Tdiff_dens();


//%------>>>>>>file : slice.c
//%
void slice_init();
void make_slice_img( int pt, double *data, long NPart, double *weight );
void field_slice( int pt, double *data, char *name, long N, double *weight );
void mag_slice();
void mach_slice();
void density_slice();
void temperature_slice();
void cren_slice();
void cree_slice();
void radio_slice();


//%------>>>>>>file : cosmology.c
//%
double E_a ( double a );
double hubble_function( double a );
double com_integ( double a, void *params );
double comoving_distance( double a );
double angular_distance( double a );
double luminosity_distance( double a );
double OmegaM( double a );
double OmegaLambda_a( double a );
double growth_factor( double a );
double growth_factor0( double a );
double PowerSpec_Efstathiou(double k);
double top_hat_filter( double k, void *params );
void top_hat_filter_k_limit( double *k0, double *k1, void *params );
double dsigma_dk( double k, void *params );
double sigma( sigma_struct ss );
void test_sigma();
double top_hat_dsigma2dmdk( double k, void *params);
double dsigma2dm( sigma_struct ss );
void MtoFilterParams( double M, void *params );
void init_ps();
double PS_dndM( double a, double M );
void test_ps();
void compute_cosmo_quantities();
void test_cos();


//%------>>>>>>file : read_params.c
//%
void read_parameters( char *fn );


//%------>>>>>>file : group_gas_ratio.c
//%
void group_gas_ratio();


