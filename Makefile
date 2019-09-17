GSL_INCL   = 
GSL_LIBS   = -lgsl -lgslcblas -lm 

HDF5_INCL  =
HDF5_OPTS  = -DH5_USE_16_API
HDF5_LIBS  = -lhdf5

FFTW_INCL  = 
FFTW_LIBS  = -ldrfftw -ldfftw

LIBS       = $(GSL_LIBS) $(HDF5_LIBS) $(FFTW_LIBS) -lm
INCL       = $(GSL_INCL) $(HDF5_INCL) $(FFTW_INCL)

OPTS       = $(HDF5_OPTS) -Wall #-O2 #-O3
DEBUG     ?=
CC         =  mpicc 

EXEC       =  ./bin/gadget-analysis
SRCS       =  $(wildcard *.c) 
MY_INCL    =  $(wildcard *.h)
OBJS       =  $(SRCS:.c=.o)

$(EXEC): $(OBJS)
	$(CC) $(OPTS) $(OBJS) $(LIBS) -o $(EXEC)

%.o:%.c $(MY_INCL) Makefile
	$(CC) $(OPTS) $(DEBUG) $(INCL) -c $< -o $@

.PHONY: clean

clean:
	-rm  $(EXEC) *.o 
	-rm  allvars.c add_params.h gadget-analysis-config.h gadget-analysis.in protos.h

test_make:
	@echo $(TEMP)
