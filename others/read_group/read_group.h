#include "stdio.h"
#include "stdlib.h"
#include "dirent.h"
#include "limits.h"
#include "string.h"

#ifdef LONGIDS
typedef unsigned long MyIDType;
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
struct group_struct *group;

