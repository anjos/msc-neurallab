# $Id: Makefile,v 1.3 2001/07/31 11:27:15 andre Exp $
# André Rabello <Andre.Rabello@ufrj.br>

# This script will install the scripts, setting their proper
# permissions in order to work with neural-lab. The configuration is
# also set.

CONFIGFILE = default.sdb
BACKUP = $(CONFIGFILE:%=%~)
BASEDIR = $(shell cd .. && pwd | sed -e "s?$(HOME)?@HOME@?" -)
TYPEDEF = "string[[:space:]]*base[[:space:]]*=[[:space:]]*"
CONFIGSTRING = $(shell cat $(CONFIGFILE) | egrep $(TYPEDEF) | \
		 sed -e s/$(TYPEDEF)//)
OTHERCONFIG = $(shell cat default.sdb | egrep @[[:print:]]*@ | grep -c file)

all: config

config: 
	@echo "[Script Makefile] Setting permissions..."
	@chmod 755 run matsh
	@echo "[Script Makefile] Configuring directories..."
	@cp $(CONFIGFILE) $(BACKUP) #make a backup
	@sed -e "s?$(CONFIGSTRING)?$(BASEDIR)?" $(BACKUP) > $(CONFIGFILE)
	@echo "[Script Makefile] Configuration done!"
	@if [ $(OTHERCONFIG) != 0 ]; then \
		echo -n "  >>> You have $(OTHERCONFIG) other config"; \
		echo " parameters that use environment variables"; \
		echo "  >>> Please, update them manually!"; \
	 fi

.PHONY: clean

GARBAGE = "*~" "*.bak" "*.pyc"
GARBDIR = .
FINDCLEANOPT = -a -type f -print0 #-maxdepth 1
clean:
	@for i in $(GARBAGE); do \
	    find $(GARBDIR) -name "$$i" $(FINDCLEANOPT) | xargs -0tr rm -f;\
	 done
