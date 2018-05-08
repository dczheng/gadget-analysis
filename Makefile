CFLAGS  =  `pkg-config --cflags gsl ` -DH5_USE_16_API
CLIBS   =  `pkg-config --libs  gsl` -lhdf5 -lm
EXEC    = ./gadget-analysis
OBJS    = main.o read_snapshot.o allvars.o analysis.o debug.o read_parameters.o\
		  read_group.o set_units.o slice.o kernel.o write_img.o\
		  auxiliary_functions.o tree.o cosmology.o ngb.o fof.o system.o\
		  conv.o
OPT     = -Wall

DEBUG ?=
CC      = mpicc $(DEBUG)

$(EXEC): $(OBJS)
	$(CC) $(OPT) $(OBJS) $(CLIBS) -o $(EXEC)
$(OBJS): allvars.h protos.h macros.h

clean: 
	rm -f $(EXEC)  $(OBJS)
