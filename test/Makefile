# Hello emacs, this is -*- makefile -*-
# $Id: Makefile,v 1.1 2001/07/13 15:48:02 andre Exp $
# André Rabello <Andre.Rabello@ufrj.br>

# This Makefile will help make in building the Calo-Neural package

# Some definition macros
CC = gcc
INCLUDES = include ../include
CFLAGS = -O2 -Wall -ansi -pedantic -D_GNU_SOURCE $(INCLUDES:%=-I%)
LIBS = m popt
LDFLAGS = -O2 -Wall $(LIBS:%=-l%)

# The name of files
SRC = src/main.c src/neural.c src/lut.c src/parameter.c ../src/real.c ../src/analysis.c ../src/event.c
OBJ = $(SRC:%.c=%.o)

all: nnet 

nnet: $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LIBS:%=-l%) -o $@

.PHONY: clean

clean:
	rm -f */*~ *~
	rm -f nnet
	rm -f $(SRC:%.c=%.o) $(SRC:%.c=%.c~) include/*~

