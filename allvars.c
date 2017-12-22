#include "allvars.h"

char sep_str[ SEP_LEN ];
long long TotNgroups, NumPart, N_Gas, SliceStart[6], SliceEnd[6];
int ThisTask, NumTask;
double BoxSize, RedShift, *KernelMat2D[6], *KernelMat3D[6];
void *CommBuffer;
char LogFile[ FILENAME_MAX ], LogBuf[500];
FILE *LogFilefd;
struct para_struct para;
struct particle_data *P;
struct sph_particle_data *SphP;

struct io_header header;
struct group_struct *group;
struct plot_struct plot_info;

#ifdef DEBUG
    float debug_f;
    int debug_i;
    long debug_l;
    double debug_d;
    char *debug_s;
#endif

gsl_integration_workspace *inte_ws;

int proj_i, proj_j;
double proj_x, proj_y, proj_size;

