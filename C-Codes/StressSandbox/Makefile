CC = gcc
CFLAGS = -g -O0 -Wall -std=c99

LIBS = -lm

INC = -I.

SRC =	../SDF-FaultRep/FuzzyFault.c \
		../FrictionLaws/Lib_NewmarkTS.c \
		../FrictionLaws/Lib_SetOfFrictionLaws.c \
		Lib_MiniVoigt.c \
		Lib_DisplacementFunctions.c \

OBJ = $(SRC:.c=.o)

HDR = $(SRC:.c=.h)

all:
	-@$(MAKE) DemoStress


DemoStress: $(OBJ) DemoStress.o
	$(CC) $(CFLAGS) -o DemoStress DemoStress.o $(INC) $(OBJ) $(LIBS)

%.o: %.c $(HDR)
	-@echo "---- Compiling $*.c -----"
	$(CC) $(CFLAGS) -c -o $*.o $*.c $(INC)


clean:
	-@rm -rf *.o
	-@rm -rf DemoStress

