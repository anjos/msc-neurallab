# $Id: Makefile,v 1.2 2001/07/31 00:39:36 andre Exp $
# André Rabello <Andre.Rabello@ufrj.br>

# This Makefile will help make in building the neural-tester package

# Some definition macros
CC = gcc
INCLUDES = include ../include
INCLUDECMD = $(INCLUDES:%=-I%)
CFLAGS = -O2 -Wall -ansi -pedantic -D_GNU_SOURCE $(INCLUDECMD)
LIBS = m popt
LDFLAGS = -O2 -Wall $(LIBS:%=-l%)
MKDEPFILE = .depend
MKDEP = makedepend -f$(MKDEPFILE)

# The name of files
SRC = src/main.c \
	src/neural.c \
	src/lut.c \
	src/parameter.c \
	../src/real.c \
	../src/analysis.c \
	../src/event.c
OBJ = $(SRC:%.c=%.o)

all: depend nnet 

nnet: $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LIBS:%=-l%) -o $@

.PHONY: clean depend

depend:
	@if [ ! -e $(MKDEPFILE) ]; then touch $(MKDEPFILE); fi 
	$(MKDEP) $(INCLUDECMD) $(SRC) $(DATASPEC:%=./src/%.c)

GARBAGE = "*~" "*.o" "TAGS" "nnet" "$(MKDEPFILE:%=%*)"
GARBDIR = . ../src $(INCLUDES)
FINDCLEANOPT = -a -type f -print0 #-maxdepth 1
clean:
	@for i in $(GARBAGE); do \
	    find $(GARBDIR) -name "$$i" $(FINDCLEANOPT) | xargs -0tr rm -f;\
	 done

tags:
	@find . $(INCLUDES) -name '*.[ch]' | etags -

# The dependencies (not obligatory)
# =================================

sinclude $(MKDEPFILE)
