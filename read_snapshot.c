#include "allvars.h"
//#define GROUP_ZDEBUG

hid_t hdf5_file, hdf5_group, hdf5_dataset, hdf5_dataspace, hdf5_dataspace_in_file,
       hdf5_dataspace_in_memory, hdf5_type, hdf5_hdf5_type_mem, hdf5_attribute, hdf5_type;
herr_t herr;
hsize_t dims[2], maxdims[2];
int ndims;
void *CommBuffer;

char *ParticleName[6] = { "Gas", "Halo", "Disk", "Bulge", "Stars", "Bndry" };

#define blockpresent_check( A )  { \
    if ( (A) ) \
        return 1; \
    else \
        return 0; \
}

int blockpresent0( enum iofields blk, int pt ) {
    switch( blk ) {
        case IO_VEL:
        case IO_POS:
        case IO_ID:
        case IO_POT:
        case IO_ACCEL:
            blockpresent_check( header.npart[pt] > 0 );

        case IO_MASS:
            blockpresent_check( (header.npart[pt] > 0 && header.mass[pt] == 0) );

        case IO_MAG:
        case IO_SFR:
        case IO_DIVB:
        case IO_DTB:
        case IO_U:
        case IO_TEMP:
        case IO_HSML:
        case IO_RHO:
        case IO_NE:
        case IO_MN:
        case IO_CR_C0:
        case IO_CR_Q0:
        case IO_CR_E0:
        case IO_CR_n0:
        case IO_CR_P0:
        case IO_CR_DISS:
        case IO_CR_DTE:
        case IO_CR_THER:
        case IO_CRE_C:
        case IO_CRE_QMIN:
        case IO_CRE_QMAX:
        case IO_CRE_ALPHA:
        case IO_CRE_N:
        case IO_CRE_E:
        case IO_DTE:
        case IO_MET:
        case IO_NHY:
        case IO_U_JUMP:
        case IO_PRE_RHO:
        case IO_PRE_U:
        case IO_PRE_XCR:
            blockpresent_check( (pt == 0 && header.npart[pt] > 0) );
        default:
            return 0;
    }
}

int blockpresent( enum iofields blk, int pt ) {
    switch( blk ) {
        case IO_POS:
        case IO_ID:
        case IO_MASS:
        case IO_RHO:
            blockpresent_check( blockpresent0( blk, pt ) );

        case IO_VEL:
#ifdef  READVEL
            blockpresent_check( blockpresent0(blk,pt));
#endif
            return 0;

        case IO_TEMP:
#ifdef READTEMP
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_HSML:
#ifdef READHSML
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_U:
#ifdef READU
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_NE:
#ifdef READELEC
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_MN:
#ifdef READMACH
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_MAG:
#ifdef READB
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_DIVB:
#ifdef READDIVB
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_SFR:
#ifdef READSFR
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_POT:
#ifdef READPOT
            blockpresent_check( blockpresent0(blk, pt) );
#endif
            return 0;

        case IO_CRE_C:
        case IO_CRE_ALPHA:
        case IO_CRE_QMIN:
        case IO_CRE_QMAX:
        case IO_CRE_N:
        case IO_CRE_E:
#ifdef READCRE
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_CR_Q0:
        case IO_CR_C0:
        case IO_CR_E0:
        case IO_CR_n0:
        case IO_CR_P0:
#ifdef READCR
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;
            
        case IO_CR_DISS:
#ifdef READCR
#ifdef READCRDISSIPATIONTIME
            blockpresent_check( blockpresent0(blk,pt) );
#endif
#endif
            return 0;
        case IO_CR_DTE:
#ifdef READCR
#ifdef READCRDTE
            blockpresent_check( blockpresent0(blk,pt) );
#endif
#endif
            return 0;

        case IO_CR_THER:
#ifdef READCR
#ifdef READCRTHERMALIZATIONTIME
            blockpresent_check( blockpresent0(blk,pt) );
#endif
#endif
            return 0;

        case IO_DTB:
#ifdef READDTB
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_ACCEL:
#ifdef READACC
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_RHO_JUMP:
#ifdef READDENSITYJUMP
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;
        case IO_MET:
#ifdef READMETALLICITY
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_NHY:
#ifdef READNEUTRALHYDROGENABUNDANCE
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;
        case IO_PRE_RHO:
        case IO_PRE_U:
        case IO_PRE_XCR:
#ifdef READPRESHOCK 
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        case IO_U_JUMP:
#ifdef READEJUMP 
            blockpresent_check( blockpresent0(blk,pt) );
#endif
            return 0;

        default:
            return 0;
    }
}

int get_block_nbytes( enum iofields blk ) {
    int block_nbytes = 0;
    switch ( blk ) {
        case IO_POS:
        case IO_VEL:
        case IO_ACCEL:
        case IO_MAG:
            block_nbytes = 3 * sizeof( OutputFloat );
            break;
        case IO_SFR:
        case IO_DIVB:
        case IO_DTB:
        case IO_MASS:
        case IO_U:
        case IO_TEMP:
        case IO_HSML:
        case IO_RHO:
        case IO_POT:
        case IO_NE:
        case IO_MN:
        case IO_MET:
        case IO_NHY:
        case IO_PRE_RHO:
        case IO_PRE_U:
        case IO_PRE_XCR:
        case IO_RHO_JUMP:
        case IO_U_JUMP:
        case IO_CR_DISS:
        case IO_CR_DTE:
        case IO_CR_THER:
        case IO_CR_C0:
        case IO_CR_Q0:
        case IO_CR_E0:
        case IO_CR_n0:
        case IO_CR_P0:
        case IO_CRE_C:
        case IO_CRE_ALPHA:
        case IO_CRE_QMIN:
        case IO_CRE_QMAX:
        case IO_CRE_N:
        case IO_CRE_E:
        case IO_DTE:
            block_nbytes = sizeof( OutputFloat );
            break;
        case IO_ID:
            block_nbytes = sizeof( MyIDType  );
            break;
    }
    return block_nbytes;
}

void get_block_dims( int pt, enum iofields blk, hsize_t (*dims)[2] ) {
    switch ( blk ) {
        case IO_POS:
        case IO_VEL:
        case IO_ACCEL:
        case IO_MAG:
            (*dims)[0] = header.npart[pt];
            (*dims)[1] = 3;
            break;
        case IO_SFR:
        case IO_DIVB:
        case IO_DTB:
        case IO_MASS:
        case IO_U:
        case IO_TEMP:
        case IO_HSML:
        case IO_RHO:
        case IO_POT:
        case IO_NE:
        case IO_ID:
        case IO_MN:
        case IO_MET:
        case IO_NHY:
        case IO_PRE_RHO:
        case IO_PRE_U:
        case IO_PRE_XCR:
        case IO_RHO_JUMP:
        case IO_U_JUMP:
        case IO_CR_DISS:
        case IO_CR_DTE:
        case IO_CR_THER:
        case IO_CR_C0:
        case IO_CR_Q0:
        case IO_CR_E0:
        case IO_CR_n0:
        case IO_CR_P0:
        case IO_CRE_C:
        case IO_CRE_ALPHA:
        case IO_CRE_QMIN:
        case IO_CRE_QMAX:
        case IO_CRE_N:
        case IO_CRE_E:
        case IO_DTE:
            (*dims)[0] = header.npart[pt];
            (*dims)[1] = 1;
            break;
    }
}

void get_dataset_name( enum iofields blk, char *buf ) {
    switch ( blk ) {
        case IO_POS:
            strcpy( buf, "Coordinates" );
            break;
        case IO_VEL:
            strcpy( buf, "Velocities" );
            break;
        case IO_ACCEL:
            strcpy( buf, "Acceleration" );
            break;
        case IO_SFR:
            strcpy( buf, "StarFormationRate" );
            break;
        case IO_MAG:
            strcpy( buf, "MagneticField" );
            break;
        case IO_DIVB:
            strcpy( buf, "DivergenceOfMagneticField" );
            break;
        case IO_DTB:
            strcpy( buf, "RateOfChangeOfMagneticField" );
            break;
        case IO_MASS:
            strcpy( buf, "Masses" );
            break;
        case IO_U:
            strcpy( buf, "InternalEnergy" );
            break;
        case IO_TEMP:
            strcpy( buf, "Temperature" );
            break;
        case IO_HSML:
            strcpy( buf, "SmoothingLength" );
            break;
        case IO_RHO:
            strcpy( buf, "Density" );
            break;
        case IO_POT:
            strcpy( buf, "Potential" );
            break;
        case IO_NE:
            strcpy( buf, "ElectronAbundance" );
            break;
        case IO_ID:
            strcpy( buf, "ParticleIDs" );
            break;
        case IO_MN:
            strcpy( buf, "MachNumber" );
            break;
        case IO_CR_C0:
            strcpy( buf, "CR_C0" );
            break;
        case IO_CR_Q0:
            strcpy( buf, "CR_q0" );
            break;
        case IO_CR_E0:
            strcpy( buf, "CR_E0" );
            break;
        case IO_CR_n0:
            strcpy( buf, "CR_n0" );
            break;
        case IO_CR_P0:
            strcpy( buf, "CR_P0" );
            break;
        case IO_CR_DISS:
            strcpy( buf, "CR_DissipationTime" );
            break;
        case IO_CR_DTE:
            strcpy( buf, "CR_DtE" );
            break;
        case IO_CR_THER:
            strcpy( buf, "CR_ThermalizationTime" );
            break;
        case IO_RHO_JUMP:
            strcpy( buf, "DensityJump" );
            break;
        case IO_U_JUMP:
            strcpy( buf, "EnergyJump" );
            break;
        case IO_MET:
            strcpy( buf, "Metallicity" );
            break;
        case IO_NHY:
            strcpy( buf, "NeutralHydrogenAbundance" );
            break;
        case IO_PRE_RHO:
            strcpy( buf, "Preshock_Density" );
            break;
        case IO_PRE_U:
            strcpy( buf, "Preshock_Energy" );
            break;
        case IO_PRE_XCR:
            strcpy( buf, "Preshock_XCR" );
            break;
        case IO_CRE_C:
            strcpy( buf, "CRE_C" );
            break;
        case IO_CRE_ALPHA:
            strcpy( buf, "CRE_Alpha" );
            break;
        case IO_CRE_QMIN:
            strcpy( buf, "CRE_qmin" );
            break;
        case IO_CRE_QMAX:
            strcpy( buf, "CRE_qmax" );
            break;
        case IO_CRE_N:
            strcpy( buf, "CRE_n" );
            break;
        case IO_CRE_E:
            strcpy( buf, "CRE_e" );
            break;
        case IO_DTE:
            strcpy( buf, "DtEnergy" );
            break;
    }
}

void get_hdf5_native_type( enum iofields blk, hid_t *hdf5_type ) {
    switch ( blk ) {
        case IO_ID:
#ifdef LONGIDS
            *hdf5_type = H5Tcopy(H5T_NATIVE_UINT64);
#else
            *hdf5_type = H5Tcopy(H5T_NATIVE_UINT);
#endif
            break;
        default:
#ifdef OUTPUT_IN_DOUBLEPRECISION
            *hdf5_type = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
            *hdf5_type = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
            break;
    }
}

void empty_buffer( enum iofields blk, int offset, int pt ) {
    int i, j;
    long n;
    OutputFloat *fp;
    MyIDType *ip;
    fp = CommBuffer;
    ip = CommBuffer;
    n = header.npart[pt];
    switch( blk ) {
        case IO_POS:
            for ( i=0; i<n; i++ ){
                for ( j=0; j<3; j++ )
                    P[offset+i].Pos[j] = *fp++;
                P[offset+i].Type = pt;
            }
            break;
        case IO_VEL:
#ifdef READVEL
            for ( i=0; i<n; i++ )
                for ( j=0; j<3; j++ )
                    P[offset+i].Vel[j] = *fp++;
#endif
            break;
        case IO_ID:
            for ( i=0; i<n; i++ )
                P[offset+i].ID = *ip++;
            break;
        case IO_MASS:
            for ( i=0; i<n; i++ )
                P[offset+i].Mass = *fp++;

            /*
            for ( i=0; i<n; i++ ) {
                if ( pt == 5 )
                    printf( "%g\n", P[offset+i].Mass );
            }
            */

            break;
        case IO_POT:
#ifdef READPOT
            for ( i=0; i<n; i++ )
                P[offset+i].Pot = *fp++;
#endif
            break;
        case IO_ACCEL:
#ifdef READACC
            for ( i=0; i<n; i++ )
                for ( j=0; j<3; j++ )
                    P[offset+i].Acc[j] = *fp++;
#endif
            break;
        case IO_RHO:
            for ( i=0; i<n; i++ )
                SphP[offset+i].Density = *fp++;
            break;
        case IO_SFR:
#ifdef READSFR
            for ( i=0; i<n; i++ )
                SphP[offset+i].sfr = *fp++;
#endif
            break;
        case IO_MAG:
#ifdef READB
            for ( i=0; i<n; i++ )
                for ( j=0; j<3; j++ )
                    SphP[offset+i].B[j] = *fp++;
#endif
            break;
        case IO_DIVB:
#ifdef READDIVB
            for ( i=0; i<n; i++ )
                SphP[offset+i].divB = *fp++;
#endif
            break;
        case IO_DTB:
#ifdef READDTB
            for ( i=0; i<n; i++ )
                SphP[offset+i].dBdt = *fp++;
#endif
            break;
        case IO_U:
#ifdef READU
            for ( i=0; i<n; i++ )
                SphP[offset+i].u = *fp++;
#endif
            break;
        case IO_TEMP:
#ifdef READTEMP
            for ( i=0; i<n; i++ )
                SphP[offset+i].Temp = *fp++;
#endif
            break;
        case IO_HSML:
#ifdef READHSML
            for ( i=0; i<n; i++ )
                SphP[offset+i].Hsml = *fp++;
#endif
            break;
        case IO_NE:
#ifdef READELEC
            for ( i=0; i<n; i++ )
                SphP[offset+i].elec = *fp++;
#endif
            break;
        case IO_MN:
#ifdef READMACH
            for ( i=0; i<n; i++ )
                SphP[offset+i].MachNumber = *fp++;
#endif
            break;
        case IO_CR_C0:
#ifdef READCR
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_C0 = *fp++;
#endif
            break;
        case IO_CR_Q0:
#ifdef READCR
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_Q0 = *fp++;
#endif
            break;
        case IO_CR_E0:
#ifdef READCR
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_E0 = *fp++;
#endif
            break;
        case IO_CR_n0:
#ifdef READCR
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_n0 = *fp++;
#endif
            break;
        case IO_CR_P0:
#ifdef READCR
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_P0 = *fp++;
#endif
        case IO_CR_DISS:
#ifdef READCR
#ifdef READCRDISSIPATIONTIME
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_DissipationTime = *fp++;
#endif
#endif
            break;
        case IO_CR_DTE:
#ifdef READCR
#ifdef READCRDTE
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_DtE = *fp++;
#endif
#endif
            break;
        case IO_CR_THER:
#ifdef READCR
#ifdef READCRTHERMALIZATIONTIME
            for ( i=0; i<n; i++ )
                SphP[offset+i].CR_ThermalizationTime = *fp++;
#endif
#endif
            break;
        case IO_RHO_JUMP:
#ifdef READDENSITYJUMP
            for ( i=0; i<n; i++ )
                SphP[offset+i].DensityJump = *fp++;
#endif
            break;
        case IO_U_JUMP:
#ifdef READENERGYJUMP
            for ( i=0; i<n; i++ )
                SphP[offset+i].EnergyJump = *fp++;
#endif
            break;
        case IO_MET:
#ifdef READMETALLICITY
            for ( i=0; i<n; i++ )
                SphP[offset+i].Metallicity = *fp++;
#endif
            break;
        case IO_NHY:
#ifdef READNEUTRALHYDROGENABUNDANCE
            for ( i=0; i<n; i++ )
                SphP[offset+i].NeutralHydrogenAbundance = *fp++;
#endif
            break;
        case IO_PRE_U:
#ifdef READPRESHOCK
            for ( i=0; i<n; i++ )
                SphP[offset+i].Preshock_Energy = *fp++;
#endif
            break;
        case IO_PRE_RHO:
#ifdef READPRESHOCK
            for ( i=0; i<n; i++ )
                SphP[offset+i].Preshock_Density = *fp++;
#endif
            break;
        case IO_PRE_XCR:
#ifdef READPRESHOCK
            for ( i=0; i<n; i++ )
                SphP[offset+i].Preshock_XCR = *fp++;
#endif
            break;
        case IO_CRE_C:
#ifdef READCRE
            for ( i=0; i<n; i++ )
                SphP[offset+i].CRE_C = *fp++;
#endif
            break;
        case IO_CRE_ALPHA:
#ifdef READCRE
            for ( i=0; i<n; i++ )
                SphP[offset+i].CRE_Alpha = *fp++;
#endif
            break;
        case IO_CRE_QMIN:
#ifdef READCRE
            for ( i=0; i<n; i++ )
                SphP[offset+i].CRE_qmin = *fp++;
#endif
            break;
        case IO_CRE_QMAX:
#ifdef READCRE
            for ( i=0; i<n; i++ )
                SphP[offset+i].CRE_qmax = *fp++;
#endif
            break;
        case IO_CRE_N:
#ifdef READCRE
            for ( i=0; i<n; i++ )
                SphP[offset+i].CRE_n = *fp++;
#endif
            break;
        case IO_CRE_E:
#ifdef READCRE
            for ( i=0; i<n; i++ )
                SphP[offset+i].CRE_e = *fp++;
#endif
            break;
        default:
            endruns( "can't occur!" );
            break;

    }
}

void read_header( char *fn ) {
    hdf5_file = H5Fopen( fn, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_group = H5Gopen(hdf5_file, "/Header");

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_ThisFile");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, header.npart);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_Total");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumPart_Total_HighWord");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "MassTable");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "BoxSize");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "NumFilesPerSnapshot");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
    H5Aclose(hdf5_attribute);

#ifndef GADGET2
        hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_IC_Info");
        H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_ic_info);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_DoublePrecision");
        H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision);
        H5Aclose(hdf5_attribute);
#endif

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Sfr");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Cooling");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_StellarAge");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Metals");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Flag_Feedback");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "HubbleParam");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Omega0");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "OmegaLambda");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
    H5Aclose(hdf5_attribute);

    hdf5_attribute = H5Aopen_name(hdf5_group, "Redshift");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
    H5Aclose(hdf5_attribute);

    H5Gclose(hdf5_group);
    H5Fclose(hdf5_file);
}

void show_header( io_header header ) {

    int i;

    put_sep0;

#define fmt "%-25s: "

    writelog( "Header Info:\n" );

#define LE_TEST
#ifdef LE_TEST
    double t;
    t = ((double)header.npart[4]) / header.npart[0];
    header.npart[0] = 1024 * 1024 * 1024 * (1-t);
    header.npart[1] = 1024 * 1024 * 1024;
    header.npart[4] = 1024 * 1024 * 1024 * t;
    header.npartTotal[0] = header.npart[0];
    header.npartTotal[1] = header.npart[1];
    header.npartTotal[4] = header.npart[4];
    header.mass[1] /= 64;
#endif

    writelog( fmt, "npart" );
    for ( i=0; i<6; i++ )
        writelog(  "%i ", header.npart[i] );
    writelog( "\n" );

    writelog( fmt, "mass" );
    for ( i=0; i<6; i++ )
        if ( header.mass[i] == 0 ) {
            writelog( "0 " );
        }
        else {
            writelog(  "%lf ", header.mass[i] );
        }
    writelog( "\n" );

    writelog( fmt, "npartTotal" );
    for ( i=0; i<6; i++ )
        writelog( "%i ", header.npartTotal[i] );
    writelog( "\n" );

    writelog( fmt, "npartTotalHighWord" );
    for ( i=0; i<6; i++ )
        writelog( "%i ", header.npartTotalHighWord[i] );
    writelog( "\n" );

#define write_header_i( a ) writelog( fmt"%i\n", &((#a)[7]), a )
#define write_header_f( a ) writelog( fmt"%lf\n", &((#a)[7]), a )
    write_header_f( header.redshift );
    write_header_f( header.HubbleParam );
    write_header_f( header.BoxSize );
    write_header_f( header.Omega0 );
    write_header_f( header.OmegaLambda );

    write_header_i( header.num_files );
    write_header_i( header.flag_sfr );
    write_header_i( header.flag_feedback );
    write_header_i( header.flag_stellarage );
    write_header_i( header.flag_metals );
    write_header_i( header.flag_entropy_instead_u );
    write_header_i( header.flag_doubleprecision );
    write_header_i( header.flag_ic_info );
    write_header_f( header.lpt_scalingfactor );
#undef write_header_i
#undef write_header_f
#undef fmt
    put_sep0;
}

void write_header( char *fn, io_header header ) {

    writelog( "write header To %s\n", fn  );
    hdf5_file = H5Fcreate( fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    //hdf5_file = H5Fopen( fn, H5F_ACC_RDWR, H5P_DEFAULT );
    hdf5_group = H5Gcreate( hdf5_file, "/Header", 0 );
    hsize_t adim[1] = { 6 };

    hdf5_dataspace = H5Screate( H5S_SIMPLE );
    H5Sset_extent_simple( hdf5_dataspace, 1, adim, NULL );
    hdf5_attribute = H5Acreate( hdf5_group, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_UINT, header.npart );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SIMPLE );
    H5Sset_extent_simple( hdf5_dataspace, 1, adim, NULL );
    hdf5_attribute = H5Acreate( hdf5_group, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SIMPLE );
    H5Sset_extent_simple( hdf5_dataspace, 1, adim, NULL );
    hdf5_attribute = H5Acreate( hdf5_group, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SIMPLE );
    H5Sset_extent_simple( hdf5_dataspace, 1, adim, NULL );
    hdf5_attribute = H5Acreate( hdf5_group, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "All.NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.num_files );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0 );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_DoublePrecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );

    hdf5_dataspace = H5Screate( H5S_SCALAR );
    hdf5_attribute = H5Acreate( hdf5_group, "Flag_IC_Info", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT );
    H5Awrite( hdf5_attribute, H5T_NATIVE_INT, &header.flag_ic_info );
    H5Aclose( hdf5_attribute );
    H5Sclose( hdf5_dataspace );
    H5Gclose( hdf5_group );
    H5Fclose( hdf5_file );
}

void free_particle_memory() {

    writelog( "free particle memory ...\n" );

#ifdef SHAREMEM
    myfree( MpiWin_P );
    myfree( MpiWin_SphP );
#else
    myfree_shared( MpiWin_P );
    myfree_shared( MpiWin_SphP );
#endif
    //MPI_Win_free( &MpiWin_P );
    //MPI_Win_free( &MpiWin_SphP );
}

void read_snapshot_test() {
    long i;

    do_sync_local( "" );
    mytimer_start();
    for( i=0; i<NumPart; i++ ) {
        if ( i % NTask_Local != ThisTask_Local )
            continue;
        P[i].Pos[0] ++;
    }
    sleep( 10 );
    do_sync_local( "" );
    mytimer_end();
    endrun( 20181212 );
}

void read_snapshot() {
    int pt, blk, nbytes, io_i;
    long i, file, offset, num;
    char file_name[MYFILENAME_MAX], buf[200], buf1[200];
    size_t BufferBytes;

    put_header( "read data" );

#ifdef OUTPUT_IN_DOUBLEPRECISION
    writelog( "****** Note: output precision is double ! ******\n" );
#else
    writelog( "****** Note: output precision is float ! ******\n" );
#endif
    mytimer_start();

    if ( ThisTask_Local == 0 ) {

        for( io_i = 0; io_i < IOGroups; io_i++ ) {
            if ( ThisTask_Master % IOGroups == io_i ) {
                //printf( "Master: %i read header ...\n", ThisTask_Master );
                //printf( "%i\n", All.StartSnapIndex );
                //endrun(0);
                sprintf( file_name, "%s_%03d.%3i.hdf5", All.FilePrefix, SnapIndex, 0 );
                if ( All.NumFilesPerSnapshot < 2 )
                    sprintf( file_name, "%s_%03d.hdf5", All.FilePrefix, SnapIndex );
                read_header( file_name );
            }
            do_sync_master( "" );
        }
    }

    MPI_Bcast( &header, sizeof( io_header ), MPI_BYTE,
            0, MpiComm_Local );

    N_Gas = NumPart = 0;
    for ( i=0; i<6; i++ ){
        NumPart += ( header.npartTotal[i] + ( ( (long) header.npartTotalHighWord[i] ) << 32 ) );
    }
    N_Gas = header.npartTotal[0] + ( ( (long)header.npartTotalHighWord[0] ) << 32 );

    for( pt=0; pt<6; pt++ ) {
        NumPart6[pt] = get_particle_num( pt );
        OffsetPart6[pt] = get_particle_offset( pt );
    }

    BoxSize = header.BoxSize;
    HalfBoxSize = BoxSize / 2;
    Redshift = header.redshift;
    Time = header.time;
    Omega0 = header.Omega0;
    OmegaLambda = header.OmegaLambda;
    OmegaBaryon = All.OmegaBaryon;
    HubbleParam = header.HubbleParam;
    writelog( "NumPart = %ld, N_Gas = %ld\n", NumPart, N_Gas );

    show_header( header );

    //printf( "ThisTask: %i, ThisTask_Local: %i, N_Gas: %li, NumPart: %li, Omega0: %g\n",
    //        ThisTask, ThisTask_Local, N_Gas, NumPart, header.Omega0 );
    //do_sync( "" );
    //endrun( 20181208 );

#ifdef SHAREMEM
    mymalloc1( P, NumPart * sizeof( ParticleData ) );
    mymalloc1( SphP, N_Gas * sizeof( SphParticleData ));
#else
    mymalloc1_shared( P, NumPart * sizeof( ParticleData ), sizeof( ParticleData ), MpiWin_P );
    mymalloc1_shared( SphP, N_Gas * sizeof( SphParticleData ), sizeof( ParticleData ), MpiWin_SphP );
#endif

    if ( ThisTask_Local == 0 ) {
        writelog( "Parallel ( %i ) read data ...\n", All.ParallelIO );
        do_sync_master( "" );
        for( io_i = 0; io_i < IOGroups; io_i++ ) {
            if ( ThisTask_Master % IOGroups == io_i ) {
                //printf( "Master: %i read data ...\n", ThisTask_Master );

                BufferBytes = header.npart[1] * sizeof( OutputFloat ) * 3;
                mymalloc1( CommBuffer, BufferBytes );

                for ( blk=0; blk<IO_NBLOCKS; blk++ ) {
                    for ( pt=0, offset=0; pt<6; pt++ ) {
                        for (file=0; file<All.NumFilesPerSnapshot; file++) {
                            if ( All.NumFilesPerSnapshot < 2 )
                                sprintf( file_name, "%s_%03d.hdf5", All.FilePrefix, SnapIndex );
                            else
                                sprintf( file_name, "%s_%03d.%3li.hdf5",
                                        All.FilePrefix, SnapIndex, file );
                                read_header( file_name );

                                if ( blockpresent( blk, pt ) ) {
                                    nbytes = get_block_nbytes( blk );
                                    if ( nbytes * header.npart[pt] > BufferBytes ){
                                        writelog( "CommBuffer IS TOO SMALL.\n" );
                                        myfree( CommBuffer );
                                        BufferBytes = nbytes * header.npart[pt];
                                        mymalloc1( CommBuffer, BufferBytes );
                                    }

                                    hdf5_file = H5Fopen( file_name, H5F_ACC_RDWR, H5P_DEFAULT );
                                    get_dataset_name( blk, buf );
                                    get_hdf5_native_type( blk, &hdf5_type );
#ifdef LE_TEST
                                    writelog( "[%5s] reading %s ...\n", ParticleName[pt], buf );
#else
                                    writelog( "[%5s] %8li reading %s ...\n", ParticleName[pt], offset, buf );
#endif
                                    sprintf( buf1, "PartType%i/%s", pt, buf );
                                    hdf5_dataset = H5Dopen( hdf5_file, buf1 );
                                    herr = H5Dread( hdf5_dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, CommBuffer );
                                    empty_buffer( blk, offset, pt );
                                    H5Dclose( hdf5_dataset );
                                    H5Tclose( hdf5_type );
                                    H5Fclose( hdf5_file );
                            }

                            offset += header.npart[pt];
                        }
                    }
                }
                myfree( CommBuffer );
            }
            do_sync_master( "" );
        }
    }

    //do_sync( "" );
    //sleep( 10000 );
    //printf( "ThisTask: %i, ThisTask_Local: %i, Position of Particle 0: (%g, %g, %g)\n",
    //        ThisTask, ThisTask_Local, P[0].Pos[0], P[0].Pos[1], P[0].Pos[2] );
    //do_sync( "" );
    //endrun( 20181209 );

    //check_data( 212121 );

    for ( pt=0, offset=0; pt<6; pt++ ) {
        num = header.npartTotal[pt] + ( ( ( long )header.npartTotalHighWord[pt] ) << 32 );
        if ( num != 0 && header.mass[pt] != 0 ){
            writelog( "set particle `%d` mass from header\n", pt );
            for ( i=offset; i<offset+num; i++ )
                P[i].Mass = header.mass[pt];
        }
        offset += num;
    }

    /*
    for ( i=0; i<N_Gas; i++ )
        for( pt=0; pt<2; pt++ )
            SphP[i].Star_BH_Num[pt] = SphP[i].Star_BH_MaxNum[pt] = 0;
            */

    //read_snapshot_test();
    mytimer_end();
#ifdef LE_TEST
    for( i=0; i<NumPart; i++ ) 
        P[i].Mass /= 64;
#endif

#ifdef READ_SNAPSHOT_TEST
    save_data_for_test();
#endif
    writelog( "read data ... done. \n" );
    //endrun( 20181210 );
}

int save_data_blockpresent( enum iofields blk, int pt ) {
    switch( blk ) {
        case IO_VEL:
        case IO_POS:
        /*
        case IO_ID:
        case IO_POT:
        case IO_ACCEL:
        */
        case IO_MASS:
        blockpresent_check( header.npart[pt] > 0 );

        case IO_MAG:
        case IO_NE:
        case IO_SFR:
        case IO_DIVB:
        case IO_U:
        case IO_TEMP:
        case IO_RHO:
        case IO_NHY:
        /*
        case IO_DTB:
        case IO_HSML:
        case IO_MN:
        case IO_CR_C0:
        case IO_CR_Q0:
        case IO_CR_E0:
        case IO_CR_n0:
        case IO_CR_P0:
        case IO_CR_DISS:
        case IO_CR_DTE:
        case IO_CR_THER:
        case IO_CRE_C:
        case IO_CRE_QMIN:
        case IO_CRE_QMAX:
        case IO_CRE_ALPHA:
        case IO_CRE_N:
        case IO_CRE_E:
        case IO_DTE:
        case IO_MET:
        case IO_U_JUMP:
        case IO_PRE_RHO:
        case IO_PRE_U:
        case IO_PRE_XCR:
        */
            blockpresent_check( (pt == 0 && header.npart[pt] > 0) );
        default:
            return 0;
    }
}

void save_data_for_test() {
    int blk, pt, i, N, offset;
    char name[100], fn[100];;
    N = 100;
    FILE *fd;
    if ( ThisTask )
        return;

    put_header( "save data for test" );

    for ( pt=0; pt<6; pt++ ) {
        if ( header.npart[pt] == 0 )
            continue;
        printf( "%s:\n", ParticleName[pt] );
        sprintf( fn, "%s.dat", ParticleName[pt] );
        fd = myfopen( "w", fn );
        offset = OffsetPart6[pt];
        for ( blk=0; blk<IO_NBLOCKS; blk++ ) {
            if ( !save_data_blockpresent( blk, pt ) )
                continue;
            get_dataset_name( blk, name );
            printf( "save %s \n", name );
            switch( blk ) {
                case IO_POS:
                    fprintf( fd, "x,y,z," );
                    break;
                case IO_VEL:
                    fprintf( fd, "Vx,Vy,vz," );
                    break;
                case IO_MAG:
                    fprintf( fd, "Bx,By,Bz," );
                    break;
                default:
                    fprintf( fd, "%s,", name );
            }
        }
        fprintf( fd, "\n" );
        for( i=offset; i<offset+N; i++ ) {
            for ( blk=0; blk<IO_NBLOCKS; blk++ ) {
                if ( ! save_data_blockpresent( blk, pt ) )
                    continue;
                get_dataset_name( blk, name );
                switch( blk ) {
                    case IO_POS:
                        fprintf( fd, "%g,%g,%g,",
                            P[i].Pos[0], 
                            P[i].Pos[1], 
                            P[i].Pos[2] 
                            );
                        break;
                    case IO_VEL:
                        fprintf( fd, "%g,%g,%g,",
                            P[i].Vel[0], 
                            P[i].Vel[1], 
                            P[i].Vel[2] 
                            );
                        break;
                    case IO_MAG:
                        fprintf( fd, "%g,%g,%g,",
                            SphP[i].B[0], 
                            SphP[i].B[1], 
                            SphP[i].B[2] 
                            );
                        break;
                    case IO_MASS:
                        fprintf( fd, "%g,", P[i].Mass );
                        break;
                    case IO_NE:
                        fprintf( fd, "%g,", SphP[i].elec );
                        break;
                    case IO_SFR:
                        fprintf( fd, "%g,", SphP[i].sfr );
                        break;
                    case IO_TEMP:
                        fprintf( fd, "%g,", SphP[i].Temp );
                        break;
                    case IO_RHO:
                        fprintf( fd, "%g,", SphP[i].Density );
                        break;
                    case IO_NHY:
                        fprintf( fd, "%g,", SphP[i].NeutralHydrogenAbundance );
                        break;
                    case IO_U:
                        fprintf( fd, "%g,", SphP[i].u );
                        break;
                }
            }
            fprintf( fd, "\n" );
        }
        fclose( fd );
    }
    put_end();
    endruns( "" );
}
