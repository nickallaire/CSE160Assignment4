SOURCES = heat2d.c heat2d_solver.c 
CC = gcc
CFLAGS = -g -Wall

default: heat2d 

heat2d_solver.o: heat2d_solver.c heat2d_solver.h 
	$(CC)  $(CFLAGS) -c heat2d_solver.c 

heat2d: heat2d_solver.o heat2d.c
	$(CC)  $(CFLAGS) -o heat2d heat2d.c heat2d_solver.o -lpthread

clean:
	-/bin/rm *o heat2d
