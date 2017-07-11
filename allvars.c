#include "allvars.h"

char para_file[ FILENAME_MAX ];
char file_prefix[ FILENAME_MAX ];
char out_file[ FILENAME_MAX ];
char out_picture_prefix[ FILENAME_MAX ];
char sep_str[ SEP_LEN ];
char group_dir[ FILENAME_MAX ];
int Num_files, box[3], TotNgroups;
float redshift, al[3], az[3], corner1[3], corner2[3], scalar_unit;
int this_task, task_num;

struct particle_struct Particle[6];
struct io_header header;
struct group_struct *group;


hid_t hdf5_file, hdf5_group, hdf5_dataset, hdf5_dataspace_in_file, hdf5_dataspace_in_memory, hdf5_dataspace, hdf5_type, hdf5_hdf5_type_mem, hdf5_attribute;
herr_t herr;
hsize_t dims[2], maxdims[2], npoints, precision;
int ndims, slice_num, slice_index_num, *slice_index, pic_xsize, pic_ysize;
