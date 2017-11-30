#include "allvars.h"

char FilePrefix[ FILENAME_MAX ];
char sep_str[ SEP_LEN ];
char GroupDir[ FILENAME_MAX ];
long long NumFiles, TotNgroups, PicSize, NumPart, N_Gas, BufferSize;
int MpcFlag;
float BoxSize, RedShift, SofteningTable[6];
void *CommBuffer;
double UnitTime_in_s,
       UnitMass_in_g,
       UnitLength_in_cm,
       UnitDensity_in_cgs,
       UnitEnergy_in_cgs,
       UnitVelocity_in_cm_per_s;
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
