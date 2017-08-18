CC      = gcc
CFLAGS  =   -DH5_USE_16_API
CLIBS   =   -lhdf5 -lm
EXEC    = ../bin/show_header
OBJS    = show_header.o
$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLIBS) -o $(EXEC)
	rm -f  $(OBJS)
%.o: %.c
	$(CC)  $(CFLAGS) -c  $<
clean: 
	rm -f $(EXEC)  $(OBJS)
