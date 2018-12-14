GSL_INCL   =  `pkg-config --cflags gsl`
GSL_LIBS   =  `pkg-config --libs  gsl`

LIBS       = $(GSL_LIBS)  -lm
INCL       = $(GSL_INCL) 

OPTS       = $(HDF5_OPTS) -Wall #-O2 #-O3
DEBUG     ?=
CC         =  gcc

EXEC       =  ./gsl_test
SRCS       =  $(wildcard *.c)
MY_INCL    =  $(wildcard *.h)
OBJS       =  $(SRCS:.c=.o)

$(EXEC): $(OBJS)
	$(CC) $(OPTS) $(OBJS) $(LIBS) -o $(EXEC)

%.o:%.c $(MY_INCL) Makefile
	$(CC) $(OPTS) $(DEBUG) $(INCL) -c $< -o $@

.PHONY: clean

clean:
	-rm  $(EXEC)  $(OBJS)
