CC      = mpicc
CFLAGS  =  `pkg-config --cflags gsl plplot ` -DH5_USE_16_API
CLIBS   =  `pkg-config --libs  gsl plplot ` -lhdf5 -lm
EXEC    = process_data
OBJS    = main.o io.o allvars.o plot.o analysis.o debug.o
$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLIBS) -o $(EXEC)
%.o: %.c
	$(CC)  $(CFLAGS) -c  $<
clean: 
	rm -f $(EXEC)  $(OBJS)
