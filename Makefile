CC = gcc
CFLAGS = -Wall -Wmissing-prototypes -O3 -std=c99

all : adensum doloops sctst

adensum : adensum.c
	$(CC) $(CFLAGS) adensum.c -o adensum.x

doloops : doloops.c
	$(CC) $(CFLAGS) doloops.c -o doloops.x

sctst : sctst.c
	$(CC) $(CFLAGS) sctst.c -o sctst.x

clean :
	rm -f *.x

