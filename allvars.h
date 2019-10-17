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
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_fit.h"
#include "mpi.h"
#include "signal.h"
#include "limits.h"
#include "sys/stat.h"
#include "macros.h"
#include "libgen.h"
#include "gadget-analysis-config.h"
#include "drfftw.h"
#include "debug.h"

#define ZDEBUG

#define MYFILENAME_MAX   FILENAME_MAX
#define MYPATHNAME_MAX   FILENAME_MAX

#define SEP_LEN  80
#define SEP_LEN0 60 
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

#define ASMTH  1.25
#define RCUT   4.5
#define NSRPTAB 1000

#define GSL_INTE_WS_LEN 10000
#define GSL_INTE_ERR_ABS  ((double)(0.0))
#define GSL_INTE_ERR_REL  ((double)(1e-3))
#define GSL_INTE_ERR_REL2 ((double)(1e-2))
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

#define GROUP_FIELD_NBLOCKS 1000
enum group_fields {
    GROUP_DENS,
    GROUP_TEMP,
    GROUP_SFR,
    GROUP_U,
    GROUP_MAG,
    GROUP_MACH,
    GROUP_CREN,
    GROUP_CREE,
    GROUP_CREC,
    GROUP_CREALPHA,
    GROUP_CREQMIN,
    GROUP_CREQMAX,
    GROUP_RAD,
    GROUP_RAD1,
    GROUP_RADINDEX
};

struct radio_inte_struct{
    double (*f) ( double, void* );
    void *params;
    double B;
    double nu;
};

#if defined(GROUPRAD) || defined(RADSLICE) || defined(TOTSPEC) || defined(GROUPSPEC) || defined(GROUPLUM) || (defined(OUTPUTGROUP) && defined(OUTPUTGROUPLUM)) || defined(RADSMOOTH)
#ifndef RADLARMOr
#define RAD
#endif
#endif

#if defined(GROUPTEMP) || defined(TEMPSLICE) || defined(PDFTDIFFDENS) || defined(PHASE) || defined(CRENTPDF) || defined(TPDF) || defined(GASRATIO) || defined(HSMLTPDf) || defined(UTPDF) || defined(BPDF) 
#ifndef COMPUTETEMP
#define COMPUTETEMP
#endif
#endif

#if defined(GROUPTEMP) || defined(GROUPU) || defined(GROUPSFR) || defined(GROUPB) || defined(GROUPMACH) || defined(GROUPCRE) || defined(GROUPRAD) || defined(GROUPSPEC) || defined(GROUPELECSPEC) || defined(GROUPTEMPPROFILE) || defined(GROUPTEMPSTACK) || defined(GROUPLUM) || defined(OUTPUTGROUP)
#ifndef FOF 
#define FOF 
#endif
#endif

#if defined(FOF) || defined(BSMOOTH) || defined(CRESMOOTH)
#ifndef TREE
#define TREE
#endif
#endif


#if defined(GROUPCRE) || defined(GROUPELECSPEC) || defined(CREPPDF) || defined(CRENSLICE) || defined(CREESLICE) || defined(CRENTPDF) || defined(RAD) || defined(RADLARMOR)
#ifndef READCRE
#define READCRE
#endif
#endif

#if defined(HSMLTPDF) || defined(HSMLDENSPDF) || defined(RADSLICE) || defined(DIVBERRDENSPDF) || defined(TOTSPEC) || defined(RAD) || defined(DIVBERRPDF) || defined(BSMOOTH) || defined(RADLARMOR)
#ifndef READHSML
#define READHSML
#endif
#endif

#if defined(RAD) || defined(GROUPB) || defined(BPDF) || defined(BDENSPDF) || defined(DIVBERRPDF) || defined(BSMOOTH) || defined(RADLARMOR)
#ifndef READB
#define READB
#endif
#endif

#if defined(GROUPSFR) || defined(BPDF) || defined(PHASE) 
#ifndef READSFR
#define READSFR
#endif
#endif

#if defined(FOF) || defined(GROUPVELDISP) || defined(GROUPKIN) 
#ifndef READVEL
#define READVEL
#endif
#endif

#if defined(GROUPU) || defined(UTPFD) || defined(GROUPCRE) || defined(CREESLICE) || defined(OUTPUTGROUPLUM) || (defined(COMPUTETEMP)&&!defined(READTEMP))
#ifndef READU
#define READU
#endif
#endif

#if defined(DIVBERRPDF) || defined(DIVBERRPDf) || defined(DIVBERRDENSPDF) 
#ifndef READDIVB
#define READDIVB
#endif
#endif

#if defined(GROUPMACH) || defined(MACHNOISE) || defined(MACHSLICE)  || defined(MACHSMOOTH)
#ifndef READMACH
#define READMACH
#endif
#endif


typedef struct ParticleData {

    double Pos[3];
    double Mass;
#ifdef READVEL
    double Vel[3];
#endif
#ifdef READACC
    double Acc[3];
#endif
    MyIDType ID;
#ifdef READPOT
    double  Pot;
#endif 
    int Type;
    short flag;
} ParticleData;

typedef struct SphParticleData {
    double Density;
#ifdef READU
    double u;
#endif
#ifdef READMACH
    double MachNumber;
#endif
#ifdef READCR
    double CR_C0;
    double CR_Q0;
    double CR_n0;
    double CR_E0;
    double CR_P0;
#endif
#ifdef READCRE
    double CRE_C;
    double CRE_Alpha;
    double CRE_qmin;
    double CRE_qmax;
    double CRE_n;
    double CRE_e;
    double CRE_P;
#endif
#ifdef READB
    double B[3];
#endif

#ifdef READDIVB
    double divB;
#endif
#ifdef READDTB
    double dBdt;
#endif
#ifdef READELEC
    double elec;
#endif
#if defined(READTEMP) || defined (COMPUTETEMP)
    double Temp;
#endif
#ifdef READHSML
    double Hsml;
#endif
#ifdef READSFR
    double sfr;
#endif

} SphParticleData;

typedef struct GlobalParams{
    char FilePrefix[ MYFILENAME_MAX ],
            FoFDir[ MYFILENAME_MAX ];

    int 
        MpcFlag,
        MFBins, BPdfBins, DivBErrPdfBins,
        GroupTempProfileRN,
        KernelInterpolation,
        ConvN,
        PowSpecNGrid, PowSpecPartType, PowSpecBins,
        ParallelIO,
        NGrid,
        DensPdfN,
        TPdfN,
        GroupTempStackRN,
        NumFilesPerSnapshot,
        GroupPotGrid,

        QNum, FreqN, FoFMinLen,
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
            SliceStartX,
            SliceStartY,
            SliceStartZ,
            SliceEndX,
            SliceEndY,
            SliceEndZ,
            UnitMass_in_g,
            UnitVelocity_in_cm_per_s,
            TreeAllocFactor, LinkLength,
            UnitLength_in_cm,
            Sigma8,
            ConvSigma, FreqMin, FreqMax, GroupMassMin,
            GroupRadFreq, GroupRadFreq1,
            OutputGroupFreq,
            RadSliceFreq,
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
            BPdfBMin, BPdfBMax,
            DivBErrPdfMin, DivBErrPdfMax,
            BDensPdfBMin, BDensPdfBMax,
            BDensPdfDensMin, BDensPdfDensMax,
            DivBerrDensPdfDivBMin, DivBerrDensPdfDivBMax,
            DivBerrDensPdfDensMin, DivBerrDensPdfDensMax,
            PosShiftX, PosShiftY, PosShiftZ, GroupSize;

} GlobalParams;

typedef struct NODE {
    double center[3];
    double len;
    long sons[8], sibling, father;
    long nextnode;
    int bitflags;
} NODE;

typedef struct group_properties{
    double mass, cm[3], vr200, vel[3], mass_table[6], size[3],
                ek, v_mean, v_disp;
    long Head, Tail, Len, npart[6];
} group_properties;

typedef struct gadget_2_cgs_unit{
    double cm, g, s, erg, density;
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
            NumPartInSlice,
            N_GasInSlice,
            *ConvKernel,
            *id_to_index, NumPart, N_Gas, OffsetPart6[6], NumPart6[6] ;
extern int 
            Proj[3], SnapIndex,
            PicSize, PicSize2,
            Ngroup, NPresentGroup,
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
            BoxSolidAngle,
            Time, Time2, Time3,
            Hubble_a, RhoBaryon,
            ComDis, AngDis, LumDis,
            Redshift, HubbleParam,
            RhoCrit, RhoM, G,
            Hubble, Omega0, OmegaLambda, OmegaBaryon,
            UnitPressure_in_cgs,
            UnitTime_in_s,
            SofteningTable[6], 
            SliceL, SliceStart[3], SliceEnd[3],
            *PartRad, *KernelMat2D[6], *KernelMat3D[6],
            *ShortRangeTablePotential;
extern gsl_integration_workspace *inte_ws;
extern gsl_integration_workspace *inte_ws2;
//extern_end

#define proj_i Proj[0]
#define proj_j Proj[1]
#define proj_k Proj[2]
#include "protos.h"
#include "signal_hander.h"
