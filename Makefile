PROG = anderson_Model2_Hurdle
OBJECTS = r1279.o nrutil.o secs.o indexx.o
MAINOBJ = $(PROG).o
MAINSRC = $(PROG).c

CC = cc -I$(HOME)/include
CFLAGS = -O

all:	$(MAINOBJ) $(OBJECTS)
	$(CC) $(MAINOBJ) $(OBJECTS) -o $(PROG) -lm

$(MAINOBJ):	$(MAINSRC) $(OBJECTS)
	$(CC) $(CFLAGS) -c $(MAINSRC)

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o
