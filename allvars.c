#include "allvars.h"

char sep_str[ SEP_LEN ];
long *id_to_index, NumPart, N_Gas;
int ThisTask, NTask;
gsl_integration_workspace *inte_ws;
FILE *LogFileFd, *MemUseFileFd;

struct NODE *Nodes, *Nodes_Base;
long MaxNodes;
long *NextNode;
long *Ngblist;

struct group_properties *Gprops;
long  *FoFNext;
int Ngroups;
char global_buffer[1024];

struct global_parameters_struct All;
struct Particle_Data *P;
struct Sph_Particle_Data *SphP;
struct io_header header;
struct image_struct image;
struct gadget_2_cgs_unit g2c;
struct malloc_struct ms;

#ifdef ZDEBUG
struct sig_struct sig;
#endif
