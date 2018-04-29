#include "allvars.h"

char sep_str[ SEP_LEN ];
long long TotNgroups, NumPart, N_Gas, SliceStart[6], SliceEnd[6], *id_to_index;
int ThisTask, NumTask;
double *KernelMat2D[6], *KernelMat3D[6], HalfBoxSize;
char LogFile[ FILENAME_MAX ], LogBuf[500];
FILE *LogFilefd;
struct global_parameters_struct All;
struct particle_data *P;
struct sph_particle_data *SphP;

struct io_header header;
struct group_struct *group;
struct image_struct image;

#ifdef DEBUG
    long debug_l[DEBUG_ARR_LEN];
    double debug_d[DEBUG_ARR_LEN];
    char debug_s[500];
#endif

gsl_integration_workspace *inte_ws;

struct NODE *Nodes, *Nodes_Base;
long MaxNodes;
long *NextNode;
long *Ngblist;

long  *fof_Next;
int Ngroups;
struct fof_info_struct *fof_info;

struct gadget_2_cgs_unit g2c;
int proj_i, proj_j, proj_k;
char malloc_var[MALLOC_VAR_NUM][MALLOC_VAR_LEN], malloc_str[100];
long long malloc_mem, malloc_var_bytes[MALLOC_VAR_NUM],
     malloc_i, malloc_n, malloc_b, malloc_max_mem;
