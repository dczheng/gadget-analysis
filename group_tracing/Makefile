CC      = gcc
CLIBS   =  -lhdf5 -lm
EXEC    = ../bin/group_tracing
OBJS    = group_tracing.o
OPT = -DLONGIDS
CFLAGS  =  -DH5_USE_16_API $(OPT)
$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLIBS) -o $(EXEC)
	rm -f  $(OBJS)
%.o: %.c
	$(CC)  $(CFLAGS) -c  $<
clean: 
	rm -f $(EXEC)  $(OBJS)
