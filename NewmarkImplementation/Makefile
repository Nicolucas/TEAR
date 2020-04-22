CC = gcc
CFLAGS = -g -O0 -Wall -std=c99

LIBS = -lm

INC = -I.

SRC =	../SDF-FaultRep/FuzzyFault.c \
		../FrictionLaws/Lib_NewmarkTS.c \
		../FrictionLaws/Lib_SetOfFrictionLaws.c \
		../StressSandbox/Lib_MiniVoigt.c \
		../StressSandbox/Lib_DisplacementFunctions.c \

OBJ = $(SRC:.c=.o)

HDR = $(SRC:.c=.h)

all:
	-@$(MAKE) StressEvolution


StressEvolution: $(OBJ) StressEvolution.o
	$(CC) $(CFLAGS) -o StressEvolution StressEvolution.o $(INC) $(OBJ) $(LIBS)

%.o: %.c $(HDR)
	-@echo "---- Compiling $*.c -----"
	$(CC) $(CFLAGS) -c -o $*.o $*.c $(INC)


clean:
	-@rm -rf *.o
	-@rm -rf StressEvolution

