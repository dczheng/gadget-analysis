#include "allvars.h"

char sep_str[ SEP_LEN ];
long long TotNgroups, *id_to_index, NumPart, N_Gas;
int ThisTask, NumTask;
gsl_integration_workspace *inte_ws;

struct NODE *Nodes, *Nodes_Base;
long MaxNodes;
long *NextNode;
long *Ngblist;

long  *fof_Next;
int Ngroups;

struct global_parameters_struct All;
struct particle_data *P;
struct sph_particle_data *SphP;
struct io_header header;
struct group_struct *group;
struct image_struct image;
struct fof_info_struct *fof_info;
struct gadget_2_cgs_unit g2c;

#ifdef DEBUG
    long debug_l[DEBUG_ARR_LEN];
    double debug_d[DEBUG_ARR_LEN];
    char debug_s[500];
#endif





char malloc_var[MALLOC_VAR_NUM][MALLOC_VAR_LEN], malloc_str[100];
long long malloc_mem, malloc_var_bytes[MALLOC_VAR_NUM],
     malloc_i, malloc_n, malloc_b, malloc_max_mem;
