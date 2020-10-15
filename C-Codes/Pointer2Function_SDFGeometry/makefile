CC = gcc
CFLAGS = -g -O0 -Wall -std=c99

LIBS = -lm

INC = -I.

SRC = Lib_SDFgeometry.c \

OBJ = $(SRC:.c=.o)

HDR = $(SRC:.c=.h)

all:
	-@$(MAKE) GetGeometryDemo

GetGeometryDemo: $(OBJ) Main_SelectFunc.o
	$(CC) $(CFLAGS) -o GetGeometryDemo Main_SelectFunc.o $(INC) $(OBJ) $(LIBS)


%.o: %.c $(HDR)
	-@echo "---- Compiling $*.c -----"
	$(CC) $(CFLAGS) -c -o $*.o $*.c $(INC)

clean:
	-@rm -rf *.o
	-@rm -rf ApplyGeometry
	-@rm -rf *.dSYM

