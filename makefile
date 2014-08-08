
SRCDIR=/home/$(USER)/NEMS/src
#installdate=latest
installdate := $(shell date '+%Y-%m-%d-%H-%M-%S')
INSTALLDIR=/home/$(USER)/ICE-INSTALLS/CICE_$(installdate)
LANLCICEGITDIR=/home/Fei.Liu/github/lanl_cice
LANLCICEDIR=/home/Fei.Liu/noscrub/lanl_cice

include $(SRCDIR)/conf/configure.nems

PWDDIR := $(shell pwd)
UTILINCS = -I$(LANLCICEDIR)/compile

MAKEFILE = makefile

LIBRARY  = libcice.a

MODULES  = cice_cap.o

MODULES_STUB  = 

DEPEND_FILES = ${MODULES:.o=.F90}

capgitname  := $(shell git remote -v | grep origin | head -1 | cut -f2 | cut -f1 -d " " )
capgithead  := $(shell git show-ref origin/master| cut -f1 -d " ")
cicegitname := $(shell cd $(LANLCICEGITDIR) && git remote -v | grep origin | head -1 | cut -f2 | cut -f1 -d " "  && cd $(PWDDIR) )
cicegithead := $(shell cd $(LANLCICEGITDIR) && git show-ref origin/master | cut -f1 -d " " && cd $(PWDDIR) )


all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULES)
	$(AR) $(ARFLAGS) $@ $?
	rm -f cice.mk.install
	@echo "# ESMF self-describing build dependency makefile fragment" > cice.mk.install
	@echo "# src location Zeus: $pwd" >> cice.mk.install
	@echo "# CICE github location:  $(cicegitname) $(cicegithead)" >> cice.mk.install
	@echo "# CICE CAP github location: $(capgitname) $(capgithead)" >> cice.mk.install
	@echo  >> cice.mk.install
	@echo "ESMF_DEP_FRONT     = cice_cap_mod" >> cice.mk.install
	@echo "ESMF_DEP_INCPATH   = $(INSTALLDIR)" >> cice.mk.install
	@echo "ESMF_DEP_CMPL_OBJS = " >> cice.mk.install
	@echo "ESMF_DEP_LINK_OBJS = $(INSTALLDIR)/libcice.a $(LANLCICEDIR)/liblanl_cice.a" >> cice.mk.install
	mkdir -p $(INSTALLDIR)
	cp -f libcice.a cice_cap_mod.mod $(INSTALLDIR) 
	cp -f cice.mk.install $(INSTALLDIR)/cice.mk

$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

$(MODULES_STUB): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

stub: $(MODULES_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_STUB)

clean:
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst depend

MKDEPENDS = $(SRCDIR)/../exe/mkDepends.pl

include $(SRCDIR)/conf/make.rules

include depend
