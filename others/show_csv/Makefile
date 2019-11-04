CC      = gcc
CFLAGS  = 
CLIBS   =  
EXEC    = ../../bin/show-csv
OBJS    = show-csv.o
$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLIBS) -o $(EXEC)
	rm -f  $(OBJS)
%.o: %.c
	$(CC)  $(CFLAGS) -c  $<
clean: 
	rm -f $(EXEC)  $(OBJS)
