GSL_INCL   = 
GSL_LIBS   = -lgsl -lgslcblas -lm 

HDF5_INCL  =
HDF5_OPTS  = -DH5_USE_16_API
HDF5_LIBS  = -lhdf5

FFTW_INCL  = 
FFTW_LIBS  = -ldrfftw -ldfftw

PYTHON    ?= $(shell which python3)

LIBS       = $(GSL_LIBS) $(HDF5_LIBS) $(FFTW_LIBS) -lm
INCL       = $(GSL_INCL) $(HDF5_INCL) $(FFTW_INCL)

OPTS       = $(HDF5_OPTS) -Wall #-O2 #-O3
DEBUG     ?=
CC         =  mpicc 

EXEC       =  ./bin/gadget-analysis
SRCS       =  $(wildcard *.c) 
MY_INCL    =  $(wildcard *.h)
OBJS       =  $(SRCS:.c=.o)

OBJS      +=  allvars.o
MY_INCL   += add_params.h protos.h

$(EXEC): $(OBJS)
	$(CC) $(OPTS) $(OBJS) $(LIBS) -o $(EXEC)

%.o:%.c $(MY_INCL) Makefile
	$(CC) $(OPTS) $(DEBUG) $(INCL) -c $< -o $@

allvars.o: allvars.h  preprocessor.py Makefile
	@echo "generate allvars.c ..."
	$(shell ${PYTHON} ./preprocessor.py 1)
	$(CC) $(OPTS) $(DEBUG) $(INCL) -c allvars.c -o allvars.o

add_params.h: allvars.h  preprocessor.py Makefile
	@echo "generate add_params.h ..."
	$(shell ${PYTHON} ./preprocessor.py 2)

protos.h: $(SRCS) preprocessor.py Makefile
	@echo "generate protos.h ..."
	$(shell ${PYTHON} ./preprocessor.py 3)

.PHONY: clean

clean:
	-rm  $(EXEC) *.o allvars.c add_params.h protos.h

test_make:
	@echo $(TEMP)
