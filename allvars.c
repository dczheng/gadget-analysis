#include "allvars.h"

char sep_str[ SEP_LEN ];
long long TotNgroups, NumPart, N_Gas;
int ThisTask, NumTask;
double BoxSize, RedShift;
void *CommBuffer;
char LogFile[ FILENAME_MAX ];
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
