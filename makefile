SRCDIR=/home/Fei.Liu/NEMS/src
include $(SRCDIR)/conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I/home/Fei.Liu/noscrub/lanl_cice/compile

LIBRARY  = libcice.a

MODULES  = cice_cap.o

MODULES_STUB  = 

DEPEND_FILES = ${MODULES:.o=.F90}

installdir := $(shell date '+%Y-%m-%d-%H-%M-%S')
capgithead := $(shell git show-ref origin/master| cut -f1 -d' ')
cicegithead := $(shell cd ../lanl_cice/ && git show-ref origin/master | cut -f1 -d' ' && cd ../lanl_cice_cap/)


all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULES)
	$(AR) $(ARFLAGS) $@ $?
	sed -e 's/timestr/$(installdir)/g' cice.mk.template > cice.mk.install && sed -i -e 's/cice_github_revision/$(cicegithead)/g' cice.mk.install && sed -i -e 's/cice_cap_github_revision/$(capgithead)/g' cice.mk.install && mkdir /home/Fei.Liu/ICE-INSTALLS/$(installdir) && cp libcice.a cice_cap_mod.mod /home/Fei.Liu/ICE-INSTALLS/$(installdir) && cp cice.mk.install /home/Fei.Liu/ICE-INSTALLS/$(installdir)/cice.mk && rm cice.mk.install
	cp libcice.a cice.mk cice_cap_mod.mod /home/Fei.Liu/cicenems
	
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
