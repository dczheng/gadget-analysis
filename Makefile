CFLAGS  =  `pkg-config --cflags gsl ` -DH5_USE_16_API
CLIBS   =  `pkg-config --libs  gsl` -lhdf5 -lm
EXEC    = ./gadget-analysis
OBJS    = main.o read_snapshot.o allvars.o analysis.o debug.o read_parameters.o\
		  set_units.o slice.o kernel.o img.o\
		  auxiliary_functions.o tree.o cosmology.o ngb.o fof.o system.o\
		  conv.o
OPT     = -Wall

DEBUG ?=
CC      = mpicc $(DEBUG) $(OPT)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLIBS) -o $(EXEC)
$(OBJS): allvars.h protos.h macros.h

clean: 
	rm -f $(EXEC)  $(OBJS)
