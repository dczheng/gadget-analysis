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
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_fit.h"
#include "mpi.h"
#include "signal.h"
#include "limits.h"
#include "sys/stat.h"
#include "macros.h"
#include "omp.h"

#define ZDEBUG

#define MYFILENAME_MAX   FILENAME_MAX
#define MYPATHNAME_MAX   FILENAME_MAX

#define SEP_LEN 80
#define GENERATIONS 8

#define SFR
#define BLACK_HOLES
#define SEC_PER_MEGAYEAR 3.155e13

#define BITFLAG_INSIDE_LINKINGLENGTH       9

#define PROTONMASS               1.6726e-24
#define ELECTRON_MASS            9.10953e-28
#define ELECTRON_CHARGE          4.8032e-10
#define BOLTZMANN                1.38066e-16
#define LIGHT_SPEED              2.9979e10
#define THOMSON_CROSS_SECTION    6.65246e-25
#define GRAVITY                  6.672e-8
#define HUBBLE                   3.2407789e-18  /* in h/sec */
#define PI                       M_PI
#define PC                       3.0857e18
#define KPC                      3.0857e21
#define MPC                      3.0857e24
#define HYDROGEN_MASSFRAC        0.76
#define GAMMA                    ( 5.0 / 3.0 )
#define GAMMA_MINUS1             ( GAMMA - 1 )

#define BCMB0                    (3.24e-6) // gauss

#define RADIO_F_INTERP

#define GSL_INTE_WS_LEN 10000
#define GSL_INTE_ERR_ABS ((double)(0.0))
#define GSL_INTE_ERR_REL ((double)(1e-4))
#define GSL_INTE_KEY GSL_INTEG_GAUSS61
extern gsl_integration_workspace *inte_ws;

#define SQR(X) ( X*X )
#define CUBE(X) ( X*X*X )

//#define LONGIDS
//#define OUTPUT_IN_DOUBLEPRECISION

#ifdef LONGIDS
typedef unsigned long MyIDType;
#else
typedef unsigned int MyIDType;
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double OutputFloat;
#else
typedef float OutputFloat;
#endif

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

#define IO_NBLOCKS 1000
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
    IO_CRE_C,
    IO_CRE_ALPHA,
    IO_CRE_QMIN,
    IO_CRE_QMAX,
    IO_CRE_N,
    IO_CRE_E,
    IO_DIVB,
    IO_DBDT,
    IO_SFR,
    IO_TEMP
};

#define GROUP_FILED_NBLOCKS 1000
enum group_fields {
    GROUP_DENS,
    GROUP_TEMP,
    GROUP_SFR,
    GROUP_MAG,
    GROUP_MACH,
    GROUP_HGEN,
    GROUP_HGEE,
    GROUP_HGEALPHA,
    GROUP_HGEQMIN,
    GROUP_HGEQMAX,
    GROUP_RAD,
    GROUP_RADP
};



typedef struct Particle_Data {
    double Pos[3];
    double Mass;
    double Vel[3];
    double Pot;
    double Acc[3];
    MyIDType ID;
    int Type;
} ParticleData;
extern ParticleData  *P;

typedef struct Sph_Particle_Data {
    double u;
    double Density;
    double Hsml;
    double NumNgb;
    double Pressure;
    double HydroAccel[3];
    double MachNumber;
    double CR_C0;
    double CR_Q0;
    double CR_n0;
    double CR_E0;
    double CRE_C;
    double CRE_Alpha;
    double CRE_qmin;
    double CRE_qmax;
    double CRE_n;
    double CRE_e;
    double *P;
    double B[3];
    double divB;
    double dBdt;
    double elec;
    double Temp;
    double BH_Mass;
    double Star_Mass;
    double sfr;
} SphParticleData;
extern SphParticleData *SphP;

extern char sep_str[ SEP_LEN ];
extern int ThisTask, NTask;
extern long *id_to_index, NumPart, N_Gas;
extern FILE *LogFileFd, *MemUseFileFd;
extern struct io_header header;
extern char global_buffer[1024];

extern struct global_parameters_struct {
    char FilePrefix[ MYFILENAME_MAX ],
         FoFDir[ MYFILENAME_MAX ],
         RadDir[ MYFILENAME_MAX ],
         GroupDir[ MYFILENAME_MAX ],
         *ToolsPath, Sproj;

    int FoF,
        ReadHge, ReadCr, ReadB, ReadMach, ReadSfr, ReadTemp,
        MpcFlag,
        Group, MF, MFBins, BPdf,
        GroupDens, GroupTemp, GroupSfr, GroupB, GroupMach, GroupHge,
        GroupRad, GroupSpec, TotSpec,
        GasState, GasDensity, GasTemperature,
        KernelInterpolation,
        ConvN, GroupEleSpec, RadSpec,
        PowSpec, PowSpecNGrid, PowSpecPartType, PowSpecBins,

        QNum, NuNum, FoFMinLen, proj_i, proj_j, proj_k,
        TreePartType, GroupIndexMin, GroupIndexMax,
        StartSnapIndex, ProjectDirection, KernelN,
        PicSize, PicSize2, NumFiles;

    double SofteningTable[6],
           UnitTime_in_s,
           UnitMass_in_g,
           UnitLength_in_cm,
           UnitDensity_in_cgs,
           UnitEnergy_in_cgs,
           UnitPressure_in_cgs,
           UnitVelocity_in_cm_per_s,
           Start[3], End[3],
           TreeAllocFactor, LinkLength,
           BoxSize, HalfBoxSize,
           *KernelMat2D[6], *KernelMat3D[6],
           Time, Time2, Time3, Hubble_a, RhoBaryon,
           ComDis, AngDis, LumDis,
           RedShift, HubbleParam, RhoCrit, RhoM, G,
           Hubble, Omega0, OmegaLambda, OmegaBaryon,
           Sigma8,
           *ConvKernel, ConvSigma, NuMin, NuMax, GroupMassMin, Freq,
           QMin, QMax, MFMmin, MFMmax, MFMSplit;

    long SliceStart[6], SliceEnd[6];
}All;

#define IMG_PROPS_START 15
#define IMG_PROPS_OTHERS ( All.PicSize - IMG_PROPS_START )
#define img_nprops     ( images.nprops[0] )
#define img_xmin       ( image.props[1] )
#define img_xmax       ( image.props[2] )
#define img_ymin       ( image.props[3] )
#define img_ymax       ( image.props[4] )
#define img_globxmin   ( image.props[5] )
#define img_globxmax   ( image.props[6] )
#define img_globymin   ( image.props[7] )
#define img_globymax   ( image.props[8] )
#define img_proj       ( image.props[9] )
#define img_z          ( image.props[10] )
#define img_min        ( image.props[11] )
#define img_max        ( image.props[12] )
#define img_globmin    ( image.props[13] )
#define img_globmax    ( image.props[14] )
#define img_props(i)   ( image.props[ IMG_PROPS_START+i ] )
struct image_struct{
    double *data, *img,
           *props;
    /* props:
     *      0 nprops
     *      1 xmin
     *      2 xmax
     *      3 ymin
     *      4 ymax
     *      5 globxmin
     *      6 globxmax
     *      7 globymin
     *      8 globymax
     *      9 projctdirection
     *      10 redshift
     *      11 min_i
     *      12 max_i
     *      13 globmin_i
     *      14 globmax_i
     *      15 ~ end others ... */

};
extern struct image_struct image;

extern struct NODE {
    double center[3];
    double len;
    long suns[8], sibling, father;
    long nextnode;
    int bitflags;
} *Nodes, *Nodes_Base;
extern long MaxNodes;
extern long *NextNode;
extern long *Ngblist;

struct group_properties{
    double mass, cm[3], vr200, vel[3], mass_table[6], size;
    long Head, Tail, Len, npart[6];
};

extern struct group_properties *Gprops;
extern long *FoFNext;
extern int Ngroups;

extern struct gadget_2_cgs_unit{
    double cm, g, s, erg;
}g2c;

extern struct aux_constants{
    double e_mec;
} aux_c;

#define MALLOC_VAR_NUM 1000
#define MALLOC_VAR_LEN 200

struct malloc_struct {
    char var[MALLOC_VAR_NUM][MALLOC_VAR_LEN], str[100];
    long mem, var_bytes[MALLOC_VAR_NUM],
       i, nn, b, max_mem;
};

extern struct malloc_struct ms;

typedef struct sigma_struct{

    double norm;
    double (*P) ( double );
    double (*filter) ( double, void* );
    void (*FilterKLimit) ( double *, double *, void * );
    void *FilterParams;

    void ( *MtoFilterParams ) ( double, void* );

    double (*dsigma2dmdk) ( double, void* );
    void *dsigma2dmdk_params;

} sigma_struct;

#include "protos.h"

#ifdef ZDEBUG
void signal_hander( int s );
void init_sig();
void empty_sig_buf();
#define ZDEBUG_NUM 10

struct sig_struct{
    char stop[500];
    char buf[ZDEBUG_NUM][500];
    int i, j;
    double d[ZDEBUG_NUM];
} sig;

#define RAISE_SIGSTOP() { \
    sprintf( sig.stop, "%s %s %i" , __FILE__, __FUNCTION__, __LINE__ ); raise( SIGSEGV );\
}

#define ZSPRINTF( k, fmt, ... ) { \
    if ( k >= ZDEBUG_NUM ) { \
        sprintf( sig.buf[0], "%d is too large.", k ); \
        RAISE_SIGSTOP(); \
    } \
    sprintf( sig.buf[k], fmt, ##__VA_ARGS__ ); \
}

#define ZCLEAR() { empty_sig_buf();  }

#else

#define ZSPRINTF( k, fmt, ... )
#define ZCLEAR()

#endif

