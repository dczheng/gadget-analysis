CC      = mpicc
CFLAGS  =  `pkg-config --cflags cfitsio gsl plplot ` -DH5_USE_16_API
CLIBS   =  `pkg-config --libs cfitsio gsl plplot ` -lhdf5 -lm
EXEC    = process_data
OBJS    = main.o io.o allvars.o plot.o analysis.o
$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLIBS) -o $(EXEC)
%.o: %.c
	$(CC)  $(CFLAGS) -c  $<
clean: 
	rm -f $(EXEC)  $(OBJS)
