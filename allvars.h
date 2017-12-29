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
#include "proto.h"

#define DEBUG

#define MAX_PARA_FILE_LINE_LEN 200
#define SEP_LEN 80
#define IO_NBLOCKS 100
//#define LONGIDS

#define GENERATIONS 8

#define SFR
#define BLACK_HOLES
#define SEC_PER_MEGAYEAR 3.155e13

#define PROTONMASS               1.6726e-24
#define ELECTRON_MASS            9.10953e-28
#define ELECTRON_CHARGE          4.8032e-10
#define BOLTZMANN                1.38066e-16
#define LIGHT_SPEED              2.9979e10
#define THOMSON_CROSS_SECTION    6.65246e-25
#define PARSEC                   3.0857e18
#define GRAVITY                  6.672e-8
#define HUBBLE                   3.24077789e-18  /* in h/sec */
#define PI M_PI

#define GSL_INTE_WS_LEN 1000
#define GSL_INTE_ERR_ABS 0.0
#define GSL_INTE_ERR_REL 1e-3
#define GSL_INTE_KEY GSL_INTEG_GAUSS61
extern gsl_integration_workspace *inte_ws;

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
    IO_CR_E0,
    IO_CR_n0,
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
    int Flag;
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
    MyFloat CR_n0;
    MyFloat CR_E0;
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
extern double *KernelMat2D[6], *KernelMat3D[6], HalfBoxSize;
extern void *CommBuffer;
extern long long NumPart, N_Gas, TotNgroups, SliceStart[6], SliceEnd[6];
extern char LogFile[ FILENAME_MAX ], LogBuf[500];
extern FILE *LogFilefd;
extern struct global_parameters_struct {
    char FilePrefix[ FILENAME_MAX ], GroupDir[ FILENAME_MAX ],
         FofFileName[ FILENAME_MAX ];
    FILE *LogFilefd;
    int StartSnapIndex, MpcFlag, ProjectDirection, KernelN,
        PicSize, BufferSize, NumFiles, HgeFlag, CrFlag, BFlag;
    double SofteningTable[6], Alpha;
    double UnitTime_in_s,
         UnitMass_in_g,
         UnitLength_in_cm,
         UnitDensity_in_cgs,
         UnitEnergy_in_cgs,
         UnitVelocity_in_cm_per_s;
    double Start[3], End[3];
    double TreeAllocFactor, LinkLength;
    double G, Hubble, Omega0, OmegaLambda, OmegaBaryon,
           BoxSize, HalfBoxSize, RedShift, HubbleParam, CriticalDensity;
    int FofMinLen;
}All;

extern struct plot_struct{
    int log_flag, istart, iend, global_colorbar_flag;
    char data_name[100],cb_label[100], xlabel[100], ylabel[100], title[100];
    double *data, h;
}plot_info;

extern struct NODE {
    double center[3];
    double len;
    long suns[8], sibling, father;
    long nextnode;
} *Nodes, *Nodes_Base;
extern long MaxNodes;
extern long *NextNode;
extern long *Ngblist;
extern struct fof_info_struct{
    long Head, Len, Tail;
    double mass, cm[3], vr200, vel[3];
} *fof_info;
extern long *fof_Next, Ngroups;

extern struct gadget_2_cgs_unit{
    double cm, g, s, erg;
}g2c;

extern int proj_i, proj_j;
extern double proj_x, proj_y, proj_size;


#ifdef DEBUG
#define DEBUG_ARR_LEN 6
    extern long debug_l[DEBUG_ARR_LEN];
    extern double debug_d[DEBUG_ARR_LEN];
    extern char debug_s[500];
#endif

