CC      = mpicc
CFLAGS  =  `pkg-config --cflags gsl` -DH5_USE_16_API
CLIBS   =  `pkg-config --libs  gsl` -lhdf5 -lm
EXEC    = g2_bin2hdf5
OBJS    = g2_bin2hdf5.o
$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLIBS) -o $(EXEC)
	rm -f $(OBJS)
%.o: %.c
	$(CC)  $(CFLAGS) -c  $<
clean: 
	rm -f $(EXEC)  $(OBJS)
