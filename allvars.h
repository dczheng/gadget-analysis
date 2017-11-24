#include "hdf5.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "unistd.h"
#include "math.h"
#include "time.h"
#include "dirent.h"
#include "gsl/gsl_integration.h"
#include "mpi.h"
#include "signal.h"
#include "unistd.h"

#define MAX_PARA_FILE_LINE_LEN 200
#define SEP_LEN 50
#define IO_NBLOCKS 100
#define LONGIDS
#define SFR
#define BLACK_HOLES
#define SEC_PER_MEGAYEAR 3.155e13
#define DEBUG

#ifdef LONGIDS
typedef unsigned long long MyIDType;
#else
typedef unsigned int MyIDType;
#endif

#ifndef DOUBLEPRECISION     /* default is single-precision */
typedef float  MyFloat;
typedef float  MyDouble;
typedef float  MyDoublePos;
#else
#if (DOUBLEPRECISION+0) == 2
typedef float   MyFloat;
typedef double  MyDouble;
typedef double  MyDoublePos;
#else
#if (DOUBLEPRECISION+0) == 3
typedef float   MyFloat;
typedef float   MyDouble;
typedef double  MyDoublePos;
#else                        /* everything double-precision */
typedef double  MyFloat;
typedef double  MyDouble;
typedef double  MyDoublePos;
#endif
#endif
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double MyOutputFloat;
#else
typedef float MyOutputFloat;
#endif

#ifdef INPUT_IN_DOUBLEPRECISION
typedef double MyInputFloat;
#else
typedef float MyInputFloat;
#endif


struct group_struct{
  int Len;
  unsigned int Offset;
  MyIDType MinID;
  MyIDType MinIDTask;
  int GrNr;
#ifndef FOF_EXTENDED_PROPERTIES
  int LenType[6];
  MyOutputFloat MassType[6];
#endif
  MyOutputFloat Mass;
  MyOutputFloat CM[3];
  MyOutputFloat Vel[3];
  MyDoublePos FirstPos[3];
#ifdef SFR
  double Sfr;
#endif
#ifdef BLACK_HOLES
  MyOutputFloat BH_Mass;
  MyOutputFloat BH_Mdot;
  MyOutputFloat MaxDens;
  int index_maxdens, task_maxdens;
#endif

#ifdef SUBFIND
  int Nsubs;
  int FirstSub;
  MyDoublePos Pos[3];
  MyOutputFloat M_TopHat200, R_TopHat200;
  MyOutputFloat M_Mean200, R_Mean200;
  MyOutputFloat M_Crit200, R_Crit200;
#ifdef SO_VEL_DISPERSIONS
  MyOutputFloat VelDisp_TopHat200, VelDisp_Mean200, VelDisp_Crit200;
#endif
#ifdef SO_BAR_INFO
  MyOutputFloat M_Mean500, R_Mean500;
  MyOutputFloat M_Crit500, R_Crit500;
  MyOutputFloat M_Crit2500, R_Crit2500;
#ifdef SO_VEL_DISPERSIONS
  MyOutputFloat VelDisp_Mean500, VelDisp_Crit500, VelDisp_Crit2500;
#endif
#endif
  int ContaminationLen;
  MyOutputFloat ContaminationMass;
#ifdef SO_BAR_INFO
  MyOutputFloat gas_mass[6], star_mass[6], temp[6], xlum[6], ygas[6];
#endif
#endif

#ifdef FOF_EXTENDED_PROPERTIES
  MyOutputFloat VelDisp, Rmax, Vmax;
  MyOutputFloat ToI[9];
  MyOutputFloat AngMom[9];
  MyOutputFloat Pos[3];
  unsigned short Origintask;
#endif

};

struct particle_struct {
    float *pos;
    float *vel;
    float *accel;
    float *mag;
    float *elec;
    MyIDType *id;
    long num;
    float *pot;
    float *m;
    float *u;
    float *rho;
    float *mn;
    float *j;
    float *c0;
    float *q0;
    float *cre_c0;
    float *cre_q0;
    float *cre_e0;
    float *cre_n0;
};

struct io_header{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  unsigned int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int flag_entropy_instead_u;
  int flag_doubleprecision;
  int flag_ic_info;
  float lpt_scalingfactor;
  char fill[18];
  char names[15][2];
};

enum iofields {
    IO_POS,
    IO_VEL,
    IO_ACCEL,
    IO_MAG,
    IO_MASS,
    IO_U,
    IO_RHO,
    IO_POT,
    IO_ELEC,
    IO_ID,
    IO_MN,
    IO_J,
    IO_CR_C0,
    IO_CR_Q0,
    IO_CRE_C0,
    IO_CRE_Q0,
    IO_CRE_E0,
    IO_CRE_n0
};


extern struct particle_struct Particle[6];
extern struct io_header header;
extern struct group_struct *group;

extern char FilePrefix[ FILENAME_MAX ];
extern char GroupDir[ FILENAME_MAX ];
extern char sep_str[ SEP_LEN ];
extern int NumFiles, TotNgroups, PicSize;
extern double UnitTime_in_s,
              UnitMass_in_g,
              UnitLength_in_cm,
              UnitDensity_in_cgs,
              UnitEnergy_in_cgs,
              UnitVelocity_in_cm_per_s,
              UnitTime_in_Megayears;

#ifdef DEBUG
    extern float debug_f;
    extern int debug_i;
    extern long debug_l;
    extern double debug_d;
    extern char debug_s[200];
#endif


extern hid_t hdf5_file, hdf5_group, hdf5_dataset, hdf5_dataspace, hdf5_dataspace_in_file, hdf5_dataspace_in_memory, hdf5_type, hdf5_hdf5_type_mem, hdf5_attribute, hdf5_type;
extern herr_t herr;
extern hsize_t dims[2], maxdims[2], npoints, precision;
extern int ndims;

void show_header( struct io_header header );
void read_header();
void read_all_data();
void free_all_memory();
void write_file( char *fn, struct io_header header, struct particle_struct *Particle);
void get_dataset_name( enum iofields blk, char *buf );
void magnetic_field_analysis();
void density_analysis();
void hg_electrons_analysis();
void group_analysis();
void read_group();
void free_group();
void signal_hander( int sig );
void endrun( int ierr );
