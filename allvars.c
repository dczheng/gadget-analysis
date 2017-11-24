#include "allvars.h"

char FilePrefix[ FILENAME_MAX ];
char sep_str[ SEP_LEN ];
char GroupDir[ FILENAME_MAX ];
int NumFiles, TotNgroups, PicSize, ndims;
double UnitTime_in_s,
       UnitMass_in_g,
       UnitLength_in_cm,
       UnitDensity_in_cgs,
       UnitEnergy_in_cgs,
       UnitVelocity_in_cm_per_s,
       UnitTime_in_Megayears;

struct particle_struct Particle[6];
struct io_header header;
struct group_struct *group;

#ifdef DEBUG
    float debug_f;
    int debug_i;
    long debug_l;
    double debug_d;
    char debug_s[200];
#endif


hid_t hdf5_file, hdf5_group, hdf5_dataset, hdf5_dataspace_in_file, hdf5_dataspace_in_memory, hdf5_dataspace, hdf5_type, hdf5_hdf5_type_mem, hdf5_attribute;
herr_t herr;
hsize_t dims[2], maxdims[2], npoints, precision;
