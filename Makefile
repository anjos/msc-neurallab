# $Id: Makefile,v 1.1 2001/07/31 00:33:52 andre Exp $
# André Rabello <Andre.Rabello@ufrj.br>

# This Makefile forms a comprehensive foundation to the use of
# neural-lab, since it can build, install and clean all distributions
# of this toolset.

# Makefile parameters
MODULE = $(shell basename `pwd`)

# Do everything
all: version train-prog test-prog script-config
	@echo "[Central Makefile] NEURAL-LAB is ready!"

# Builds the default test program
test-prog:
	@echo "[Central Makefile] Building test program"
	$(MAKE) -C test

# Builds the default train program
train-prog:
	@echo "[Central Makefile] Building train program"
	$(MAKE) -C train


.PHONY = version clean clean-train clean-test clean-script dist shot tags

# Check permissions on scripting directory and change some
# configuration parameters
script-config:
	@echo "[Central Makefile] Installing scripts"
	$(MAKE) -C script

# The version rule
version:
	@echo \* This is NEURAL-LAB version `cat VERSION`;
	@echo \* Andre Rabello \<Andre\.Rabello\@mail\.com\>;

dist: clean
	@echo \* Creating distribution. Current date will be written on DATE...
	@echo \* Today is `date +%A,\ %d\ of\ %B\ of\ %Y`
	@echo `date` > DATE
	@cd ..; \
	 tar cvf - $(MODULE) | gzip > $(MODULE)-`cat $(MODULE)/VERSION`.tar.gz

shot: clean
	@echo \* Creating a snapshot of today\'s source...
	@echo \* The date will be written on DATE.
	@echo \* Today is `date +%A,\ %d\ of\ %B\ of\ %Y`
	@echo `date` > DATE
	@cd ..; tar cvf - $(MODULE) | gzip > $(MODULE)-`date +%Y.%m.%d`.tar.gz

tags:
	@echo \* Creating [x]emacs TAGS file...
	@find . -name '*.[ch]' | etags -

# Things to be cleaned-up
clean-train:
	$(MAKE) -C train clean

clean-test:
	$(MAKE) -C test clean

clean-script:
	$(MAKE) -C script clean

GARBAGE = "*~" "DATE"
GARBDIR = .
FINDCLEANOPT = -a -type f -maxdepth 1 -print0
clean: clean-train clean-test clean-script
	@for i in $(GARBAGE); do \
	    find $(GARBDIR) -name "$$i" $(FINDCLEANOPT) | xargs -0tr rm -f; \
	 done
