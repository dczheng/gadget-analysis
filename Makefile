CC      = mpicc
CFLAGS  =  `pkg-config --cflags gsl plplot ` -DH5_USE_16_API
CLIBS   =  `pkg-config --libs  gsl plplot ` -lhdf5 -lm
EXEC    = ./bin/gadget-analysis
OBJS    = main.o io.o  allvars.o plot.o analysis.o debug.o
$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLIBS) -o $(EXEC)
	rm -f  $(OBJS)
%.o: %.c
	$(CC)  $(CFLAGS) -c  $<
clean: 
	rm -f $(EXEC)  $(OBJS)
