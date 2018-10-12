
CFLAGS  :=  `pkg-config --cflags gsl`  -DH5_USE_16_API
CLIBS   :=  `pkg-config --libs  gsl`   -lhdf5 -lm
OPT     :=  -Wall # -O3
CC      :=  mpicc
DEBUG ?=

EXEC    :=  ./gadget-analysis

SRCS    :=  $(wildcard *.c)
HS      :=  $(wildcard *.h)
OBJS    :=  $(SRCS:.c=.o)

$(EXEC): $(OBJS)
	$(CC) $(DEBUG) $(OPT) $(OBJS) $(CLIBS) -o $(EXEC)

$(OBJS): $(HS)

.PHONY: clean

clean:
	-rm  $(EXEC)  $(OBJS)
