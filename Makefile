GSL_INCL   =  `pkg-config --cflags gsl`
GSL_LIBS   =  `pkg-config --libs  gsl`

HDF5_INCL  =
HDF5_OPTS  = -DH5_USE_16_API
HDF5_LIBS  = -lhdf5

FFTW_INCL  = -I/mnt/ddnfs/data_users/dczheng/local/usr/fftw-2.1.5_mpich/include
FFTW_LIBS  = -L/mnt/ddnfs/data_users/dczheng/local/usr/fftw-2.1.5_mpich/lib -ldrfftw -ldfftw

PYTHON ?= $(shell which python3)

LIBS       = $(GSL_LIBS) $(HDF5_LIBS) $(FFTW_LIBS) -lm
INCL       = $(GSL_INCL) $(HDF5_INCL) $(FFTW_INCL)

OPTS       = $(HDF5_OPTS) -Wall #-O2 #-O3
DEBUG     ?=
CC         =  mpicc 

EXEC       =  ./bin/gadget-analysis
SRCS       =   allvars.c analysis.c check_flag.c conv.c correlation.c cosmology.c cre.c debug.c field.c fof.c gas_ratio.c grid.c group_analysis.c img.c kernel.c main.c mf.c mymath.c ngb.c part_info.c part_radio.c pdf.c phase.c powerspec.c pre_proc.c radio.c read_params.c read_snapshot.c set_units.c slice.c system.c temp.c total_radio_spectrum.c tree.c
MY_INCL    =  allvars.h macros.h protos.h add_params.h
OBJS       =  $(SRCS:.c=.o)

$(EXEC): $(OBJS)
	$(CC) $(OPTS) $(OBJS) $(LIBS) -o $(EXEC)

%.o:%.c $(MY_INCL) Makefile
	$(CC) $(OPTS) $(DEBUG) $(INCL) -c $< -o $@

allvars.c: allvars.h preprocessor.py
	$(shell ${PYTHON} ./preprocessor.py)

add_params.h: allvars.h preprocessor.py
protos.h: $(SRCS) preprocessor.py


.PHONY: clean

clean:
	-rm  $(EXEC)  $(OBJS) allvars.c add_params.h protos.h
