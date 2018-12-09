GSL_INCL   =  `pkg-config --cflags gsl`
GSL_LIBS   =  `pkg-config --libs  gsl`

FFTW_INCL  = -I/mnt/ddnfs/data_users/dczheng/local/usr/fftw-2.1.5_mpich/include
FFTW_LIBS  = -L/mnt/ddnfs/data_users/dczheng/local/usr/fftw-2.1.5_mpich/lib -lsrfftw -lsfftw

LIBS       = $(HDF5_LIBS) $(FFTW_LIBS) -lm
INCL       = $(HDF5_INCL) $(FFTW_INCL)

OPTS       = $(HDF5_OPTS) -Wall -fopenmp  #-O3

HYBRID     ?= ""

ifeq ($(HYBRID), "")
CC         = gcc
else
CC         = mpicc
OPTS      += -DHYBRID
endif

EXEC       =  ./hybrid_test
SRCS       =  $(wildcard *.c)
MY_INCL    =  $(wildcard *.h)
OBJS       =  $(SRCS:.c=.o)

$(EXEC): $(OBJS)
	$(CC) $(OPTS) $(OBJS) $(LIBS) -o $(EXEC)

%.o:%.c $(MY_INCL)
	$(CC) $(OPTS) $(DEBUG) $(INCL) -c $< -o $@

.PHONY: clean

clean:
	-rm  $(EXEC)  $(OBJS)
