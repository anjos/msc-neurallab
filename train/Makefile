# Hello emacs, this is -*- makefile -*-
# $Id: Makefile,v 1.1 2001/07/13 15:48:02 andre Exp $
# André Rabello <Andre.Rabello@ufrj.br>

# This Makefile will help make in building the Calo-Neural package

# Some definition macros
CC = gcc
INCLUDES = include ../include
CFLAGS = -O3 -Wall -ansi -pedantic -D_GNU_SOURCE $(INCLUDES:%=-I%)
LIBS = f2c m popt
LDFLAGS = -O3 -Wall $(LIBS:%=-l%)

# The name of files
SRC = src/jetnet34.c ../src/event.c ../src/analysis.c src/parameter.c src/main.c 
OBJ = $(SRC:%.c=%.o)

all: train

train: $(OBJ)
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: clean

clean:
	rm -f */*~ *~
	rm -f train
	rm -f include/*~ ../include/*~ $(SRC:%.c=%.o) $(SRC:%.c=%.c~)

