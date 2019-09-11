#include "hdf5.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "unistd.h"

#define SEP_LEN 100

#define NUMCRPOP 1

struct io_header2 {
  int npart[6];                        /*!< number of particles of each type in this file */
  double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
                                            stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                         /*!< time of snapshot file */
  double redshift;                     /*!< redshift of snapshot file */
  int flag_sfr;                        /*!< flags whether the simulation was including star formation */
  int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];          /*!< total number of particles of each type in this snapshot. This can be
                                            different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;                    /*!< flags whether cooling was included  */
  int num_files;                       /*!< number of files in multi-file snapshot */
  double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                       /*!< matter density in units of critical density */
  double OmegaLambda;                  /*!< cosmological constant parameter */
  double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
  int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
  unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  char fill[60];	               /*!< fills to 256 Bytes */
};

struct io_header3_no_cr{
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

struct io_header3_cr{
  int npart[6];
  double mass[6];
  double SpectralIndex_CR_Pop[ NUMCRPOP ];
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
  char fill[18-8*NUMCRPOP];
  char names[15][2];
};

struct io_header3_no_cr header3_no_cr;
struct io_header3_cr header3_cr;
struct io_header2 header2;
char sep_str[ SEP_LEN ];

void read_header3_no_cr( char *fn ) {

    //fprintf( stdout, "Reading header3_no_cr From %s\n", File_Name );
    hid_t hdf5_file, hdf5_group, hdf5_attribute;
    hdf5_file = H5Fopen( fn, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_group = H5Gopen(hdf5_file, "/Header");

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_ThisFile");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, header3_no_cr.npart);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_Total");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header3_no_cr.npartTotal);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_Total_HighWord");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header3_no_cr.npartTotalHighWord);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "MassTable");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header3_no_cr.mass);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_no_cr.time);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "BoxSize");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_no_cr.BoxSize);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumFilesPerSnapshot");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_no_cr.num_files);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_IC_Info");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_no_cr.flag_ic_info);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_DoublePrecision");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_no_cr.flag_doubleprecision);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Cooling");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_no_cr.flag_cooling);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Feedback");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_no_cr.flag_feedback);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Metals");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_no_cr.flag_metals);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Sfr");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_no_cr.flag_sfr);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_StellarAge");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_no_cr.flag_stellarage);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "HubbleParam");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_no_cr.HubbleParam);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Omega0");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_no_cr.Omega0);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "OmegaLambda");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_no_cr.OmegaLambda);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_no_cr.time);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Redshift");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_no_cr.redshift);
    H5Aclose(hdf5_attribute);

    H5Gclose(hdf5_group);
    H5Fclose(hdf5_file);
}

void read_header3_cr( char *fn ) {

    //fprintf( stdout, "Reading header3_cr From %s\n", File_Name );
    hid_t hdf5_file, hdf5_group, hdf5_attribute;
    hdf5_file = H5Fopen( fn, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_group = H5Gopen(hdf5_file, "/Header");

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_ThisFile");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, header3_cr.npart);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_Total");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header3_cr.npartTotal);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_Total_HighWord");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header3_cr.npartTotalHighWord);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "MassTable");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header3_cr.mass);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_cr.time);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "BoxSize");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_cr.BoxSize);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumFilesPerSnapshot");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_cr.num_files);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_IC_Info");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_cr.flag_ic_info);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_DoublePrecision");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_cr.flag_doubleprecision);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Cooling");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_cr.flag_cooling);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Feedback");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_cr.flag_feedback);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Metals");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_cr.flag_metals);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Sfr");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_cr.flag_sfr);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_StellarAge");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header3_cr.flag_stellarage);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "HubbleParam");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_cr.HubbleParam);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Omega0");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_cr.Omega0);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "OmegaLambda");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_cr.OmegaLambda);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_cr.time);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Redshift");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header3_cr.redshift);
    H5Aclose(hdf5_attribute);

    H5Gclose(hdf5_group);
    H5Fclose(hdf5_file);
}

void read_header2( char *fn ) {

    //fprintf( stdout, "Reading header2 From %s\n", File_Name );
    hid_t hdf5_file, hdf5_group, hdf5_attribute;
    hdf5_file = H5Fopen( fn, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_group = H5Gopen(hdf5_file, "/Header");

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_ThisFile");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, header2.npart);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_Total");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header2.npartTotal);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_Total_HighWord");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header2.npartTotalHighWord);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "MassTable");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header2.mass);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header2.time);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "BoxSize");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header2.BoxSize);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumFilesPerSnapshot");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header2.num_files);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Cooling");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header2.flag_cooling);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Feedback");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header2.flag_feedback);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Metals");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header2.flag_metals);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Sfr");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header2.flag_sfr);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_StellarAge");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header2.flag_stellarage);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "HubbleParam");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header2.HubbleParam);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Omega0");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header2.Omega0);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "OmegaLambda");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header2.OmegaLambda);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header2.time);
    H5Aclose(hdf5_attribute);

    H5Gclose(hdf5_group);
    H5Fclose(hdf5_file);
}

void show_header2( struct io_header2 header ) {
    int i;
    fputs( sep_str, stdout );
    fputs( "Header2 Info: \n", stdout );
    fprintf( stdout, "%-25s: ", "npart" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npart[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "mass" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%lf ", header.mass[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "npartTotal" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npartTotal[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "npartTotalHighWord" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npartTotalHighWord[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: %lf\n", "redshift",             header.redshift );
    fprintf( stdout, "%-25s: %lf\n", "time",                  header.time );
    fprintf( stdout, "%-25s: %i\n", "flag_sfr",               header.flag_sfr );
    fprintf( stdout, "%-25s: %i\n", "flag_feedback",          header.flag_feedback );
    fprintf( stdout, "%-25s: %i\n", "num_files",              header.num_files );
    fprintf( stdout, "%-25s: %lf\n", "BoxSize",               header.BoxSize );
    fprintf( stdout, "%-25s: %lf\n", "Omega0",                header.Omega0 );
    fprintf( stdout, "%-25s: %lf\n", "OmegaLambda",           header.OmegaLambda );
    fprintf( stdout, "%-25s: %lf\n", "Hubbleparam",           header.HubbleParam );
    fprintf( stdout, "%-25s: %i\n", "flag_stellarge",         header.flag_stellarage );
    fprintf( stdout, "%-25s: %i\n", "flag_metals",            header.flag_metals );
    fprintf( stdout, "%-25s: %i\n", "flag_entropy_instead_u", header.flag_entropy_instead_u );
    fputs( sep_str, stdout );
}

void show_header3_cr( struct io_header3_cr header ) {
    int i;
    fputs( sep_str, stdout );
    fputs( "Header3 Info: \n", stdout );
    fprintf( stdout, "%-25s: ", "npart" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npart[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "mass" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%lf ", header.mass[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "npartTotal" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npartTotal[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "npartTotalHighWord" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npartTotalHighWord[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: %lf\n", "redshift",             header.redshift );
    fprintf( stdout, "%-25s: %lf\n", "time",                  header.time );
    fprintf( stdout, "%-25s: %i\n", "flag_sfr",               header.flag_sfr );
    fprintf( stdout, "%-25s: %i\n", "flag_feedback",          header.flag_feedback );
    fprintf( stdout, "%-25s: %i\n", "num_files",              header.num_files );
    fprintf( stdout, "%-25s: %lf\n", "BoxSize",               header.BoxSize );
    fprintf( stdout, "%-25s: %lf\n", "Omega0",                header.Omega0 );
    fprintf( stdout, "%-25s: %lf\n", "OmegaLambda",           header.OmegaLambda );
    fprintf( stdout, "%-25s: %lf\n", "Hubbleparam",           header.HubbleParam );
    fprintf( stdout, "%-25s: %i\n", "flag_stellarge",         header.flag_stellarage );
    fprintf( stdout, "%-25s: %i\n", "flag_metals",            header.flag_metals );
    fprintf( stdout, "%-25s: %i\n", "flag_entropy_instead_u", header.flag_entropy_instead_u );
    fprintf( stdout, "%-25s: %i\n", "flag_doubleprecision",   header.flag_doubleprecision );
    fprintf( stdout, "%-25s: %i\n", "flag_ic_info",           header.flag_ic_info );
    fprintf( stdout, "%-25s: %f\n", "lpt_scalingfactor",      header.lpt_scalingfactor );
    fprintf( stdout, "%-25s: %lf\n", "SpectralIndex_CR_Pop[0]",      header.SpectralIndex_CR_Pop[0] );
    fputs( sep_str, stdout );
}

void show_header3_no_cr( struct io_header3_no_cr header ) {
    int i;
    fputs( sep_str, stdout );
    fputs( "Header3 Info: \n", stdout );
    fprintf( stdout, "%-25s: ", "npart" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npart[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "mass" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%lf ", header.mass[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "npartTotal" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npartTotal[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: ", "npartTotalHighWord" );
    for ( i=0; i<6; i++ )
        fprintf( stdout, "%i ", header.npartTotalHighWord[i] );
    fprintf( stdout, "\n" );
    fprintf( stdout, "%-25s: %lf\n", "redshift",             header.redshift );
    fprintf( stdout, "%-25s: %lf\n", "time",                  header.time );
    fprintf( stdout, "%-25s: %i\n", "flag_sfr",               header.flag_sfr );
    fprintf( stdout, "%-25s: %i\n", "flag_feedback",          header.flag_feedback );
    fprintf( stdout, "%-25s: %i\n", "num_files",              header.num_files );
    fprintf( stdout, "%-25s: %lf\n", "BoxSize",               header.BoxSize );
    fprintf( stdout, "%-25s: %lf\n", "Omega0",                header.Omega0 );
    fprintf( stdout, "%-25s: %lf\n", "OmegaLambda",           header.OmegaLambda );
    fprintf( stdout, "%-25s: %lf\n", "Hubbleparam",           header.HubbleParam );
    fprintf( stdout, "%-25s: %i\n", "flag_stellarge",         header.flag_stellarage );
    fprintf( stdout, "%-25s: %i\n", "flag_metals",            header.flag_metals );
    fprintf( stdout, "%-25s: %i\n", "flag_entropy_instead_u", header.flag_entropy_instead_u );
    fprintf( stdout, "%-25s: %i\n", "flag_doubleprecision",   header.flag_doubleprecision );
    fprintf( stdout, "%-25s: %i\n", "flag_ic_info",           header.flag_ic_info );
    fprintf( stdout, "%-25s: %f\n", "lpt_scalingfactor",      header.lpt_scalingfactor );
    fputs( sep_str, stdout );
}

int main( int argc, char** argv ) {
    if ( argc < 3 ) {
        fprintf( stderr, "Too few command line parameters! \n" );
        return 0;
    }
    int hdf5_flag, g_flag, cr_flag, dummy;
    FILE *fd;
    memset( sep_str, '-', SEP_LEN-1 );
    sep_str[ SEP_LEN-1 ] = '\n';
    hdf5_flag = 1;
    if ( NULL == strstr( argv[1], "hdf5" ) ) hdf5_flag = 0;
    g_flag = atoi( argv[2] );
    cr_flag = 0;
    if ( argc == 4 ) cr_flag = 1;
    fprintf( stdout, "hdf5_flag = %i\n", hdf5_flag );
    fprintf( stdout, "g_flag = %i\n", g_flag );
    fprintf( stdout, "cr_flag = %i\n", cr_flag );
    switch ( g_flag ) {
        case 2:
            if ( hdf5_flag ) {
                read_header2( argv[1] );
                show_header2( header2 );
            }
            else {
                fd = fopen( argv[1], "r" );
                fread( &dummy, sizeof( dummy ), 1, fd );
                fread( &header2, sizeof( header2 ), 1, fd );
                fread( &dummy, sizeof( dummy ), 1, fd );
                fclose( fd );
                show_header2( header2 );
            }
            break;
        case 3:
            if ( hdf5_flag ) {
                if ( cr_flag ) {
                    read_header3_cr( argv[1] );
                    show_header3_cr( header3_cr );
                }
                else{
                    read_header3_no_cr( argv[1] );
                    show_header3_no_cr( header3_no_cr );
                }
            }else{
                fd = fopen( argv[1], "r" );
                fread( &dummy, sizeof( dummy ), 1, fd );
                if ( cr_flag ) {
                    fread( &header3_cr, sizeof( header3_cr ), 1, fd );
                    show_header3_cr( header3_cr );
                }
                else{
                    fread( &header3_no_cr, sizeof( header3_no_cr ), 1, fd );
                    show_header3_no_cr( header3_no_cr );
                }
                fread( &dummy, sizeof( dummy ), 1, fd );
                fclose( fd );
            }
            break;
        default:
            fprintf( stderr, "g_flag must be set 2 or 3 !\n" );
    }
}

