# Objects
# OBJ =

CC = gcc
SCR = ag.c
CFLAGS = -g -Wall
CLIBS = -lm

ag: ag.c ag.h
	$(CC) $(CFLAGS) $(CLIBS) ag.c -o ag

clean:
	rm -rf ag
