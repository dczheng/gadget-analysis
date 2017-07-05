CC      = gcc
CFLAGS  =   -DH5_USE_16_API
CLIBS   =   -lhdf5 -lm
EXEC    = show_header
OBJS    = show_header.o
$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLIBS) -o $(EXEC)
%.o: %.c
	$(CC)  $(CFLAGS) -c  $<
clean: 
	rm -f $(EXEC)  $(OBJS)
