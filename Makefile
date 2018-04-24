CC      = mpicc
CFLAGS  =  `pkg-config --cflags gsl ` -DH5_USE_16_API
CLIBS   =  `pkg-config --libs  gsl giza` -lhdf5 -lm
EXEC    = ./gadget-analysis
OBJS    = main.o read_snapshot.o allvars.o analysis.o debug.o read_parameters.o\
		  read_group.o set_units.o slice.o projection.o kernel.o plot.o\
		  auxiliary_functions.o tree.o cosmology.o ngb.o fof.o
OPT     = -Wall

$(EXEC): $(OBJS)
	$(CC) $(OPT) $(OBJS) $(CLIBS) -o $(EXEC)
	rm -f  $(OBJS)
$(OBJS): allvars.h

clean: 
	rm -f $(EXEC)  $(OBJS)
