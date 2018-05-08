#include "allvars.h"

char sep_str[ SEP_LEN ];
long TotNgroups, *id_to_index, NumPart, N_Gas;
int ThisTask, NumTask;
gsl_integration_workspace *inte_ws;
FILE *LogFileFd;

struct NODE *Nodes, *Nodes_Base;
long MaxNodes;
long *NextNode;
long *Ngblist;

struct fof_properties *FoFProps;
long  *FoFNext;
int Ngroups;

struct global_parameters_struct All;
struct particle_data *P;
struct sph_particle_data *SphP;
struct io_header header;
struct group_struct *group;
struct image_struct image;
struct fof_info_struct *fof_info;
struct gadget_2_cgs_unit g2c;
struct malloc_struct ms;

#ifdef ZDEBUG
struct sig_struct sig;
#endif
