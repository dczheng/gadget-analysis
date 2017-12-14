#include "hdf5.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "unistd.h"
#include "math.h"
#include "time.h"
#include "dirent.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf_gamma.h"
#include "mpi.h"
#include "signal.h"
#include "unistd.h"
#include "giza.h"
#include "limits.h"
#include "sys/stat.h"

#define DEBUG

#define MAX_PARA_FILE_LINE_LEN 200
#define SEP_LEN 80
#define IO_NBLOCKS 100
#define LONGIDS
#define SFR
#define BLACK_HOLES
#define SEC_PER_MEGAYEAR 3.155e13

#define PROTONMASS   1.6726e-24
#define ELECTRONMASS 9.10953e-28
#define ELECTRONCHARGE  4.8032e-10
#define BOLTZMANN      1.38066e-16

#define cgs_cm ( header.HubbleParam / para.UnitLength_in_cm )
#define cgs_g  ( header.HubbleParam / para.UnitMass_in_g )
#define cgs_s  ( header.HubbleParam / para.UnitTime_in_s )
#define cgs_erg ( cgs_g*cgs_cm*cgs_cm/(cgs_s*cgs_s))
#define cgs_keV (1.602e-9*cgs_erg)
#define deg 1.0
#define cgs_mp (PROTONMASS * cgs_g)
#define cgs_me (ELECTRONMASS * cgs_g)
#define k_B (BOLTZMANN * cgs_erg / deg)
#define cgs_c (2.9979e10*cgs_cm/cgs_s)
#define c 2.9979e10
#define HBAR ( 1.05457e-27 * cgs_cm * cgs_cm * cgs_g / cgs_s )
#define e2  ( HBAR * c / 137.04 )
#define statcoul sqrt( cgs_erg * cgs_cm )
#define cgs_mpc2 ( cgs_mp * cgs_c * cgs_c )
#define cgs_mec2 ( cgs_me * cgs_c * cgs_c )
#define cgs_c2   ( c * c )

#define GSL_INTE_WS_LEN 1000
#define GSL_INTE_ERR_ABS 0.0
#define GSL_INTE_ERR_REL 1e-3
#define GSL_INTE_KEY GSL_INTEG_GAUSS15

#define SQR(X) ( X*X )
#define CUBE(X) ( X*X*X )

#ifdef LONGIDS
typedef unsigned long long MyIDType;
#else
typedef unsigned int MyIDType;
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double MyOutputFloat;
typedef double MyFloat;
typedef double MyDouble;
typedef double MyDoublePos;
#else
typedef float MyOutputFloat;
typedef float MyFloat;
typedef float MyDouble;
typedef float MyDoublePos;
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
    IO_ID,
    IO_POS,
    IO_MASS,
    IO_VEL,
    IO_ACCEL,
    IO_MAG,
    IO_NE,
    IO_U,
    IO_RHO,
    IO_POT,
    IO_MN,
    IO_CR_C0,
    IO_CR_Q0,
    IO_CRE_C0,
    IO_CRE_Q0,
    IO_CRE_E0,
    IO_CRE_n0,
    IO_DIVB,
    IO_DBDT
};


extern struct io_header header;
extern struct group_struct *group;

extern struct particle_data {
    MyFloat Pos[3];
    MyFloat Mass;
    MyFloat Vel[3];
    MyFloat Pot;
    MyFloat Acc[3];
    MyIDType ID;
    int Type;
} *P;

extern struct sph_particle_data {
    MyFloat Entropy;
    MyFloat Density;
    MyFloat Hsml;
    MyFloat NumNgb;
    MyFloat Pressure;
    MyFloat HydroAccel[3];
    MyFloat MachNumber;
    MyFloat CR_C0;
    MyFloat CR_Q0;
    MyFloat CRE_C0;
    MyFloat CRE_Q0;
    MyFloat CRE_E0;
    MyFloat CRE_n0;
    MyFloat B[3];
    MyFloat divB;
    MyFloat dBdt;
    MyFloat elec;
    double vL;
    double P;
} *SphP;

extern char sep_str[ SEP_LEN ];
extern int ThisTask, NumTask;
extern double BoxSize, RedShift, *KernelMat2D[6], *KernelMat3D[6];
extern void *CommBuffer;
extern long long NumPart, N_Gas, TotNgroups, SliceStart[6], SliceEnd[6];
extern char LogFile[ FILENAME_MAX ], LogBuf[200];
extern FILE *LogFilefd;
extern struct para_struct {
    char FilePrefix[ FILENAME_MAX ], GroupDir[ FILENAME_MAX ];
    FILE *LogFilefd;
    int StartSnapIndex, MpcFlag, ProjectDirection, KernelN,
        PicSize, BufferSize, NumFiles;
    double SofteningTable[6], Alpha;
    double UnitTime_in_s,
         UnitMass_in_g,
         UnitLength_in_cm,
         UnitDensity_in_cgs,
         UnitEnergy_in_cgs,
         UnitVelocity_in_cm_per_s;
    double StartX, EndX, StartY, EndY, StartZ, EndZ;
}para;

extern int cpn;
extern double *cp, *red, *green, *blue, affine[6];
extern char cb_s, cb_label[100], title[100], xlabel[100], ylabel[100];
extern gsl_integration_workspace *inte_ws;
extern FILE *fp_tmp;
extern int proj_i, proj_j;
extern double proj_x, proj_y, proj_size;


#ifdef DEBUG
    extern float debug_f;
    extern int debug_i;
    extern long debug_l;
    extern double debug_d;
    extern char debug_s[200];
#endif


void show_header( struct io_header header );
void read_header();
void read_snapshot();
void free_memory();
void get_dataset_name( enum iofields blk, char *buf );
void magnetic_field_analysis();
void group_analysis();
void read_group();
void free_group();
void signal_hander( int sig );
void endrun( int ierr );
void read_parameters();
void set_units();
void slice();
void init_projection();
void print_log( char *log );

void analysis();
void pos_analysis();
void mach_analysis();
void gas_density_analysis();
void hg_electrons_analysis();
void cosmic_rays_analysis();
void magnetic_field_analysis();
void radio_radiation_analysis();
double kernel( double q );
void init_kernel_matrix();
void free_kernel_matrix();
