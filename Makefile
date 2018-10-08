CFLAGS  =  `pkg-config --cflags gsl`  -DH5_USE_16_API
CLIBS   =  `pkg-config --libs  gsl`   -lhdf5 -lm
EXEC    =  ./gadget-analysis
OBJS    =  main.o read_snapshot.o allvars.o analysis.o debug.o read_parameters.o\
		   set_units.o slice.o kernel.o img.o\
		   auxiliary_functions.o tree.o cosmology.o ngb.o fof.o system.o\
		   conv.o gas_state.o gas_temp.o total_radio_spectrum.o\
		   mf.o part_radio.o gas_density.o group_analysis.o\
		   mymath.o
OPT     =  -Wall # -O3

DEBUG ?=
CC      =  mpicc $(DEBUG) $(OPT)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLIBS) -o $(EXEC)
$(OBJS): allvars.h protos.h macros.h

clean:
	rm -f $(EXEC)  $(OBJS)
