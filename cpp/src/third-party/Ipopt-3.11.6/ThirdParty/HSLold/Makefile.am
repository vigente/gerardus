# Copyright (C) 2003, 2009 International Business Machines and others
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 2312 2011-06-19 14:29:40Z stefan $

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = INSTALL.HSL

if BUILD_COINHSL

lib_LTLIBRARIES = libcoinhsl.la

libcoinhsl_la_SOURCES =

if COIN_HAS_MA27
  libcoinhsl_la_SOURCES += ma27ad.f
endif

if COIN_HAS_MA57
  libcoinhsl_la_SOURCES += ma57ad.f
endif

if COIN_HAS_MC19
  libcoinhsl_la_SOURCES += mc19ad.f
endif

if COIN_HAS_MA28
  libcoinhsl_la_SOURCES += ma28ad.f ma30ad.f mc13dd.f mc20ad.f \
	mc21ad.f mc22ad.f mc23ad.f mc24ad.f
endif

if COIN_HAS_MA77
  libcoinhsl_la_SOURCES += hsl_ma77d.f90 hsl_ma77d.deps90.f90 hsl_ma77d.deps77.f
endif

libcoinhsl_la_LDFLAGS = $(LT_LDFLAGS)

# dependencies
hsl_ma77d.lo: hsl_ma77d.deps90.lo

#############################################################################
#  Generating a shared library for delay loading

if EXTRA_SHARED_COMMAND

$(HSL_SHARED_LIB): $(thirdpartylib_LTLIBRARIES) $(DEFFILE)
	bla=; \
	for i in $(nodist_libcoinhsl_la_SOURCES); do \
	  case $$i in \
	    *.f) \
	      obj=`echo $$i| sed -e 's|\(.*\).[cf]|.libs/\1.$(OBJEXT)|g'` ;\
	      if test -r $$obj; then \
	        bla="$$bla $$obj" ;\
	      else \
	        bla="$$bla `echo $$i| sed -e 's|\(.*\).[cf]|\1.$(OBJEXT)|g'`" ;\
	      fi ;\
	      ;; \
	  esac ;\
	done ;\
	$(F77_SHARED_CMD) $$bla -o $(HSL_SHARED_LIB) $(F77_DLL_LINKFLAGS)

if NEED_DEFFILE
$(DEFFILE): config.status
	echo 'LIBRARY LIBHSL' > $(DEFFILE)
	echo 'EXPORTS' >> $(DEFFILE)
	for i in $(EXPORT_SYMBOLS); do \
	  echo "  $$i" >> $(DEFFILE); \
	done
endif

endif  # EXTRA_SHARED_COMMAND

includecoindir = $(includedir)/coin/ThirdParty
install-exec-local: $(HSL_SHARED_LIB)
if EXTRA_SHARED_COMMAND
	test -z "$(DESTDIR)$(libdir)/coin/ThirdParty" || $(mkdir_p) "$(DESTDIR)$(libdir)/coin/ThirdParty"
	$(libLTLIBRARIES_INSTALL) $(HSL_SHARED_LIB) "$(DESTDIR)$(libdir)/coin/ThirdParty/$(HSL_SHARED_LIB)"
endif
	$(install_sh_DATA) config_coinhsl.h $(DESTDIR)$(includecoindir)/CoinHslConfig.h

uninstall-local:
	rm -f "$(DESTDIR)$(libdir)/coin/ThirdParty/$(HSL_SHARED_LIB)"
	rm -f $(DESTDIR)$(includecoindir)/CoinHslConfig.h

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = coinhsl.pc

endif  # BUILD_COINHSL

test:
	@echo "No test available for HSL."

# Need to overwrite automake definition since otherwise libtool does not
# understand the fortran file (need to add --tag=F77)
LTFCCOMPILE = $(LIBTOOL) --tag=F77 --mode=compile $(FC) $(AM_FCFLAGS) $(FCFLAGS)


########################################################################
#                         Maintainer Stuff                             #
########################################################################

CLEANFILES = $(HSL_SHARED_LIB) $(DEFFILE)

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
