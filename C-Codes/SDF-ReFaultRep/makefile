
CC = gcc
CFLAGS = -g -O0 -Wall -std=c99

LIBS = -lm

INC = -I.

SRC =	FuzzyFault.c \
	sdf.c \

OBJ = $(SRC:.c=.o)

HDR = $(SRC:.c=.h)

all:
	-@$(MAKE) fuzzyfaultdemo1
	-@$(MAKE) sdfdemo1

fuzzyfaultdemo1: $(OBJ) demo1.o
	$(CC) $(CFLAGS) -o fuzzyfaultdemo1 demo1.o $(INC) $(OBJ) $(LIBS)

sdfdemo1: $(OBJ) demo2.o
	$(CC) $(CFLAGS) -o sdfdemo1 demo2.o $(INC) $(OBJ) $(LIBS)

%.o: %.c $(HDR)
	-@echo "---- Compiling $*.c -----"
	$(CC) $(CFLAGS) -c -o $*.o $*.c $(INC)

clean:
	-@rm -rf *.o
	-@rm -rf fuzzyfaultdemo1
	-@rm -rf sdfdemo1
	-@rm -rf *.dSYM

