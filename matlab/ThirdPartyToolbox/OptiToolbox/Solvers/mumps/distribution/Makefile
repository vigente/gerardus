#
#  This file is part of MUMPS 4.10.0, built on Tue May 10 12:56:32 UTC 2011
#
topdir = .
libdir = $(topdir)/lib

default:	dexamples

.PHONY: default alllib all s d c z \
	sexamples dexamples cexamples zexamples \
	mumps_lib requiredobj libseqneeded clean

alllib:		s d c z
all:		sexamples dexamples cexamples zexamples

s:
	$(MAKE) ARITH=s mumps_lib
d:
	$(MAKE) ARITH=d mumps_lib
c:
	$(MAKE) ARITH=c mumps_lib
z:
	$(MAKE) ARITH=z mumps_lib


# Is Makefile.inc available ?
Makefile.inc:
	@echo "######################################################################"
	@echo "# BEFORE COMPILING MUMPS, YOU SHOULD HAVE AN APPROPRIATE FILE"
	@echo "# Makefile.inc AVALAIBLE. PLEASE LOOK IN THE DIRECTORY ./Make.inc FOR"
	@echo "# EXAMPLES OF Makefile.inc FILES, AT Make.inc/Makefile.inc.generic"
	@echo "# IN CASE YOU NEED TO BUILD A NEW ONE AND READ THE MAIN README FILE"
	@echo "######################################################################"
	@exit 1

include Makefile.inc

mumps_lib: requiredobj
	(cd src ; $(MAKE) $(ARITH))

sexamples:	s
	(cd examples ; $(MAKE) s)

dexamples:	d
	(cd examples ; $(MAKE) d)

cexamples:	c
	(cd examples ; $(MAKE) c)

zexamples:	z
	(cd examples ; $(MAKE) z)


requiredobj: Makefile.inc $(LIBSEQNEEDED) $(libdir)/libpord$(PLAT)$(LIBEXT)

# dummy MPI library (sequential version)

libseqneeded:
	(cd libseq; $(MAKE))

# Build the libpord.a library and copy it into $(topdir)/lib
$(libdir)/libpord$(PLAT)$(LIBEXT):
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cd $(LPORDDIR); \
	  $(MAKE) CC="$(CC)" CFLAGS="$(OPTC)" AR="$(AR)" RANLIB="$(RANLIB)" OUTC=$(OUTC) LIBEXT=$(LIBEXT); \
	fi;
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cp $(LPORDDIR)/libpord$(LIBEXT) $@; \
	fi;

clean:
	(cd src; $(MAKE) clean)
	(cd examples; $(MAKE) clean)
	(cd $(libdir); $(RM) *$(PLAT)$(LIBEXT))
	(cd libseq; $(MAKE) clean)
	if [ $(LPORDDIR) != "" ] ; then \
	  cd $(LPORDDIR); $(MAKE) realclean; \
        fi;

