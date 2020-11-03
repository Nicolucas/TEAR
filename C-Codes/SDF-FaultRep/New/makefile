
CC = gcc
CFLAGS = -g -O0 -Wall -std=c99

LIBS = -lm

INC = -I.

SRC =	 sdf_FaultGeom.c \

OBJ = $(SRC:.c=.o)

HDR = $(SRC:.c=.h)

all:
	-@$(MAKE) FaultGeom

FaultGeom: $(OBJ) Demo.o
	$(CC) $(CFLAGS) -o FaultGeom Demo.o $(INC) $(OBJ) $(LIBS)



%.o: %.c $(HDR)
	-@echo "---- Compiling $*.c -----"
	$(CC) $(CFLAGS) -c -o $*.o $*.c $(INC)

clean:
	-@rm -rf *.o
	-@rm -rf FaultGeom
	-@rm -rf *.dSYM

