CFLAGS  =  `pkg-config --cflags gsl ` -DH5_USE_16_API
CLIBS   =  `pkg-config --libs  gsl` -lhdf5 -lm
EXEC    = ./gadget-analysis
OBJS    = main.o read_snapshot.o allvars.o analysis.o debug.o read_parameters.o\
		  read_group.o set_units.o slice.o kernel.o io.o\
		  auxiliary_functions.o tree.o cosmology.o ngb.o fof.o system.o
OPT     = -Wall

DEBUG ?=
CC      = mpicc $(DEBUG)

$(EXEC): $(OBJS)
	$(CC) $(OPT) $(OBJS) $(CLIBS) -o $(EXEC)
$(OBJS): allvars.h

clean: 
	rm -f $(EXEC)  $(OBJS)
