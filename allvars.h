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
#include "libgen.h"

#define ZDEBUG

#define MYFILENAME_MAX   FILENAME_MAX
#define MYPATHNAME_MAX   FILENAME_MAX

#define SEP_LEN 80
#define SEP_LEN0 30
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

#define MEC2                     8.187103009502731e-07
#define LS2                      8.98740441e+20


#define BCMB0                    (3.24e-6) // gauss

#define GSL_INTE_WS_LEN 10000
#define GSL_INTE_ERR_ABS ((double)(0.0))
#define GSL_INTE_ERR_REL ((double)(1e-4))
#define GSL_INTE_KEY GSL_INTEG_GAUSS61

#define SQR(X) ( (X)*(X) )
#define CUBE(X) ( SQR(X) * (X) )

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

typedef struct io_header{
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
} io_header;

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
    IO_CR_P0,
    IO_CRE_C,
    IO_CRE_ALPHA,
    IO_CRE_QMIN,
    IO_CRE_QMAX,
    IO_CRE_N,
    IO_CRE_E,
    IO_DIVB,
    IO_DBDT,
    IO_SFR,
    IO_HSML,
    IO_TEMP
};

#define GROUP_FILED_NBLOCKS 1000
enum group_fields {
    GROUP_DENS,
    GROUP_TEMP,
    GROUP_SFR,
    GROUP_MAG,
    GROUP_MACH,
    GROUP_CREN,
    GROUP_CREE,
    GROUP_CREC,
    GROUP_CREALPHA,
    GROUP_CREQMIN,
    GROUP_CREQMAX,
    GROUP_RAD,
    GROUP_RADP
};

struct radio_inte_struct{
    double (*f) ( double, void* );
    void *params;
    double B;
    double nu;
};


typedef struct ParticleData {
    double Pos[3];
    double Mass;
    double Vel[3];
    double Pot;
    double Acc[3];
    MyIDType ID;
    int Type;
} ParticleData;

typedef struct SphParticleData {
    double Density;
    double u;

    /*
    double Hsml;
    double NumNgb;
    double HydroAccel[3];
    */

    double MachNumber;
    double CR_C0;
    double CR_Q0;
    double CR_n0;
    double CR_E0;
    double CR_P0;
    double CRE_C;
    double CRE_Alpha;
    double CRE_qmin;
    double CRE_qmax;
    double CRE_n;
    double CRE_e;
    double CRE_P;
    //double *Rad;
    double B[3];
    double divB;
    double dBdt;
    double elec;
    double Temp;
    double Hsml;
    /*
    int Star_BH_Num[2];
    int Star_BH_MaxNum[2];
    MyIDType *Star_BH_Index[2];
    */
    double sfr;
} SphParticleData;

typedef struct GlobalParams{
    char FilePrefix[ MYFILENAME_MAX ],
            FoFDir[ MYFILENAME_MAX ],
         RadDir[ MYFILENAME_MAX ];

    int FoF,
        ReadCre, ReadCr, ReadB, ReadMach, ReadSfr, ReadTemp, ReadVel,
        ReadElec, Readu, ReadHsml,
        MpcFlag,
        Group, MF, MFBins, BPdf,
        GroupTemp, GroupSfr, GroupB, GroupMach, GroupCre,
        GroupTempProfileRN,
        GroupTempProfile,
        MachSlice,
        BSlice, UnitAreaSlice, CREnSlice, RadSlice,
        GroupRad, GroupSpec, TotSpec,
        Phase, DensitySlice, TemperatureSlice,
        KernelInterpolation,
        ConvN, GroupEleSpec, RadSpec,
        PowSpec, PowSpecNGrid, PowSpecPartType, PowSpecBins,
        CrePressurePdf, TabF, Tree, ParallelIO,
        CorrTdiffDens,
        PdfTdiffDens,
        NGrid,
        CorrGas, CorrDM,
        DensPdf, DensPdfN,
        TPdf, TPdfN,
        CrenTPdf,
        GroupTempStack,
        GroupTempStackRN,
        GroupFixedSize,
        GroupGasRatio,
        GasRatio,
        FieldCrenTDens,
        HsmlTPdf,
        UTPdf,
        HsmlDensPdf,
        Gadget2,
        NumFilesPerSnapshot,

        QNum, NuNum, FoFMinLen,
        TreePartType,
        StartSnapIndex, ProjectDirection, KernelN,
        PicSize;

    double 
            SofteningGas   ,
            SofteningHalo  ,
            SofteningDisk  ,
            SofteningBulge ,
            SofteningStar  ,
            SofteningBndry ,
            StartX,
            StartY,
            StartZ,
            EndX,
            EndY,
            EndZ,
            UnitMass_in_g,
            UnitVelocity_in_cm_per_s,
            TreeAllocFactor, LinkLength,
            UnitLength_in_cm,
            Sigma8,
            ConvSigma, NuMin, NuMax, GroupMassMin, Freq,
            QMin, QMax, MFMmin, MFMmax, MFMSplit,
            PhaseTempMin, 
            PhaseTempMax,
            PhaseDensMax,
            PhaseDensMin,
            DensPdfMin, DensPdfMax,
            OmegaBaryon,
            TPdfMin, TPdfMax,
            GroupTempStackRmin,
            GroupTempStackRmax,
            GroupTempProfileRmin,
            GroupTempProfileRmax,
            CrenTPdfTMin,
            CrenTPdfTMax,
            CrenTPdfnMin,
            CrenTPdfnMax,
            FieldCrenTDensDensMin,
            FieldCrenTDensDensMax,
            FieldCrenTDensTMin,
            FieldCrenTDensTMax,
            PosShiftX, PosShiftY, PosShiftZ, GroupSize;

} GlobalParams;

typedef struct NODE {
    double center[3];
    double len;
    long suns[8], sibling, father;
    long nextnode;
    int bitflags;
} NODE;

typedef struct group_properties{
    double mass, cm[3], vr200, vel[3], mass_table[6], size[3];
    long Head, Tail, Len, npart[6];
} group_properties;

typedef struct gadget_2_cgs_unit{
    double cm, g, s, erg;
} gadget_2_cgs_unit;

typedef struct physical_constants_in_gadget_unit{
    double e, m_e, m_p, c, G, mec2, c2, e_mec, sigma_t;
}physical_constants_in_gadget_unit;

typedef struct physical_constants_in_cgs_unit{
    double e, m_e, m_p, c, G, mec2, c2, e_mec, sigma_t;
} physical_constants_in_cgs_unit;

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



#define IMG_PROPS_START 8 
#define IMG_PROPS_OTHERS ( All.PicSize - IMG_PROPS_START )
#define img_z          ( image.props[0] )
#define img_xmin       ( image.props[1] )
#define img_xmax       ( image.props[2] )
#define img_xlog       ( image.props[3] )
#define img_ymin       ( image.props[4] )
#define img_ymax       ( image.props[5] )
#define img_ylog       ( image.props[6] )
#define img_proj       ( image.props[7] )
#define img_props(i)   ( image.props[ IMG_PROPS_START+i ] )
typedef struct image_struct{
    double *img,  *num, *img_tmp, *num_tmp,
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

}image_struct;

#define MALLOC_VAR_NUM 1000
#define MALLOC_VAR_LEN 200

typedef struct malloc_struct {
    char var[MALLOC_VAR_NUM][MALLOC_VAR_LEN], str[100], mem_str[100];
    long mem, var_bytes[MALLOC_VAR_NUM],
       i, nn, b, max_mem;
} malloc_struct;


//extern_start
extern SphParticleData *SphP;
extern ParticleData  *P;
extern GlobalParams All;
extern io_header header;
extern image_struct image;
extern gadget_2_cgs_unit g2c;
extern physical_constants_in_cgs_unit guc;
extern physical_constants_in_cgs_unit cuc;
extern NODE *Nodes, *Nodes_Base;
extern image_struct image;
extern group_properties *Gprops;
extern malloc_struct ms;

extern long 
            MaxNodes, *NextNode, *Ngblist,
            *FoFNext,
            SliceStart[6], SliceEnd[6],
            *ConvKernel,
            *id_to_index, NumPart, N_Gas, OffsetPart6[6], NumPart6[6] ;
extern int 
            proj_i, proj_j, proj_k, SnapIndex,
            PicSize, PicSize2,
            Ngroups,
            ThisTask, NTask, MasterTask, ThisTask_Local, NTask_Local,
            ThisTask_Master, NTask_Master, IOGroups,
            WinDisp,
            NumThreadsPerSnapshot;

extern char 
            Sproj,
         GroupDir[ MYFILENAME_MAX ],
            OutputDir[ MYFILENAME_MAX ],
            sep_str[ SEP_LEN ], sep_str0[ SEP_LEN0 ],
            *ToolsPath ;

extern MPI_Win 
            MpiWin_P, MpiWin_SphP, MpiWin_PartRad ;
extern MPI_Comm  
            MpiComm_Local, MpiComm_Master ;
extern MPI_Aint WinSize;
extern FILE 
            *LogFileFd, *UsedMemFileFd ;
extern double 
            UnitDensity_in_cgs,
            UnitEnergy_in_cgs,
            BoxSize, HalfBoxSize,
            Time, Time2, Time3,
            Hubble_a, RhoBaryon,
            ComDis, AngDis, LumDis,
            Redshift, HubbleParam,
            RhoCrit, RhoM, G,
            Hubble, Omega0, OmegaLambda,
            UnitPressure_in_cgs,
            UnitTime_in_s,
            SofteningTable[6], Start[6], End[6],
            *PartRad, *KernelMat2D[6], *KernelMat3D[6] ;
extern gsl_integration_workspace *inte_ws;
//extern_end

#include "protos.h"
#include "auxfuns.h"

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


