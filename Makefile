CC      = mpicc
CFLAGS  =  `pkg-config --cflags gsl ` -DH5_USE_16_API
CLIBS   =  `pkg-config --libs  gsl giza` -lhdf5 -lm
EXEC    = ./gadget-analysis
OBJS    = main.o read_snapshot.o  allvars.o analysis.o debug.o read_parameters.o\
		  read_group.o set_units.o slice.o print_log.o\
		  group_analysis.o\
		  mach_analysis.o\
		  gas_density_analysis.o\
		  pos_analysis.o\
		  magnetic_field_analysis.o\
		  cosmic_ray_analysis.o\
		  hge_electrons_analysis.o\
		  radio_radiation_analysis.o

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLIBS) -o $(EXEC)
	rm -f  $(OBJS)
%.o: %.c
	$(CC)  $(CFLAGS) -c  $<
clean: 
	rm -f $(EXEC)  $(OBJS)
