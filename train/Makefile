# $Id: Makefile,v 1.2 2001/07/31 00:39:41 andre Exp $
# André Rabello <Andre.Rabello@ufrj.br>

# This Makefile will help make in building the Calo-Neural package

# Some definition macros
CC = gcc
INCLUDES = include ../include
INCLUDECMD = $(INCLUDES:%=-I%)
CFLAGS = -O3 -Wall -ansi -pedantic -D_GNU_SOURCE $(INCLUDES:%=-I%)
LIBS = f2c m popt
LDFLAGS = -O3 -Wall $(LIBS:%=-l%)
MKDEPFILE = .depend
MKDEP = makedepend -f$(MKDEPFILE)

# The name of files
SRC = src/jetnet34.c \
	../src/event.c \
	../src/analysis.c \
	src/parameter.c \
	src/main.c 
OBJ = $(SRC:%.c=%.o)

all: depend train

depend:
	@if [ ! -e $(MKDEPFILE) ]; then touch $(MKDEPFILE); fi 
	$(MKDEP) $(INCLUDECMD) $(SRC) $(DATASPEC:%=./src/%.c)

train: $(OBJ)
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: clean depend

GARBAGE = "*~" "*.o" "TAGS" "train" "$(MKDEPFILE:%=%*)"
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
