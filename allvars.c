#include "allvars.h"

char sep_str[ SEP_LEN ];
long long TotNgroups, NumPart, N_Gas, SliceStart[6], SliceEnd[6];
int ThisTask, NumTask;
double BoxSize, RedShift;
void *CommBuffer;
char LogFile[ FILENAME_MAX ], LogBuf[200];
FILE *LogFilefd;
struct para_struct para;
struct particle_data *P;
struct sph_particle_data *SphP;

struct io_header header;
struct group_struct *group;

#ifdef DEBUG
    float debug_f;
    int debug_i;
    long debug_l;
    double debug_d;
    char debug_s[200];
#endif

int cpn;
double *cp, *red, *green, *blue, affine[6];
char cb_s, cb_label[100], title[100], xlabel[100], ylabel[100];
gsl_integration_workspace *inte_ws;
FILE *fp_tmp;

