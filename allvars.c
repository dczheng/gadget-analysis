#include "allvars.h"

char sep_str[ SEP_LEN ], sep_str0[ SEP_LEN0 ];
long *id_to_index, NumPart, N_Gas, NumPart6[6], OffsetPart6[6];
int ThisTask, NTask, MasterTask, ThisTask_Local, NTask_Local, ThisTask_Master, NTask_Master,
    IOGroups;
MPI_Comm  MpiComm_Local, MpiComm_Master;
MPI_Win MpiWin_P, MpiWin_SphP, MpiWin_PartRad;
MPI_Aint WinSize;
int WinDisp;

double *PartRad;

gsl_integration_workspace *inte_ws;
FILE *LogFileFd, *UsedMemFileFd;

struct NODE *Nodes, *Nodes_Base;
long MaxNodes;
long *NextNode;
long *Ngblist;

struct group_properties *Gprops;
long  *FoFNext;
int Ngroups;

struct global_parameters_struct All;
ParticleData *P;
SphParticleData *SphP;
struct io_header header;
struct image_struct image;
struct gadget_2_cgs_unit g2c;
struct physical_constants_in_gadget_unit guc;
struct physical_constants_in_cgs_unit cuc;
struct malloc_struct ms;

#ifdef ZDEBUG
struct sig_struct sig;
#endif
