GSL_INCL   =  `pkg-config --cflags gsl`
GSL_LIBS   =  `pkg-config --libs  gsl`

HDF5_INCL  =
HDF5_OPTS  = -DH5_USE_16_API
HDF5_LIBS  = -lhdf5

FFTW_INCL  = -I/mnt/ddnfs/data_users/dczheng/local/usr/fftw-2.1.5_mpich/include
FFTW_LIBS  = -L/mnt/ddnfs/data_users/dczheng/local/usr/fftw-2.1.5_mpich/lib -ldrfftw -ldfftw

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
	-rm  $(EXEC)  $(OBJS)
