# Copyright (C) 2003, 2009 International Business Machines and others
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 2696 2012-12-24 14:00:33Z stefan $

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = INSTALL.HSL

if BUILD_COINHSL

includecoindir = $(includedir)/coin/ThirdParty
includecoin_HEADERS =

lib_LTLIBRARIES = libcoinhsl.la

libcoinhsl_la_SOURCES = coinhsl/common/deps.f

if COIN_HAS_MC19
  libcoinhsl_la_SOURCES += coinhsl/mc19/mc19d.f
endif

if COIN_HAS_MA27
  libcoinhsl_la_SOURCES += coinhsl/ma27/ma27d.f
endif

if COIN_HAS_MA28
  libcoinhsl_la_SOURCES += coinhsl/ma28/ma28d.f
endif

if COIN_HAS_MA57
  libcoinhsl_la_SOURCES += coinhsl/ma57/ma57d.f
endif

if COIN_HAS_DEPSF90
  libcoinhsl_la_SOURCES += coinhsl/common/deps90.f90
endif

if COIN_HAS_HSL_MC68
  libcoinhsl_la_SOURCES += coinhsl/hsl_mc68/C/hsl_mc68i_ciface.f90
  includecoin_HEADERS += coinhsl/hsl_mc68/C/hsl_mc68i.h
endif

if COIN_HAS_HSL_MA77
  libcoinhsl_la_SOURCES += coinhsl/hsl_ma77/hsl_ma77d.f90 coinhsl/hsl_ma77/C/hsl_ma77d_ciface.f90
  includecoin_HEADERS += coinhsl/hsl_ma77/C/hsl_ma77d.h
endif

if COIN_HAS_HSL_MA86
  libcoinhsl_la_SOURCES += coinhsl/hsl_ma86/hsl_ma86d.f90 coinhsl/hsl_ma86/C/hsl_ma86d_ciface.f90
  includecoin_HEADERS += coinhsl/hsl_ma86/C/hsl_ma86d.h
endif

if COIN_HAS_HSL_MA97
  libcoinhsl_la_SOURCES += coinhsl/hsl_ma97/hsl_ma97d.f90 coinhsl/hsl_ma97/C/hsl_ma97d_ciface.f90
  includecoin_HEADERS += coinhsl/hsl_ma97/C/hsl_ma97d.h
endif

if COIN_HAS_METIS
else
  libcoinhsl_la_SOURCES += metis_dummy.c
endif

libcoinhsl_la_LDFLAGS = $(LT_LDFLAGS)

# dependencies
#  need to define dependencies in terms of variables, since automake does not allow multiple rules for the same target
ma77=hsl_ma77d
ma86=hsl_ma86d
ma97=hsl_ma97d
mc68=hsl_mc68i
$(ma77).lo: deps90.lo
$(ma86).lo: deps90.lo
$(ma97).lo: deps90.lo
$(ma77)_ciface.lo: $(ma77).lo
$(ma86)_ciface.lo: $(ma86).lo
$(ma97)_ciface.lo: $(ma97).lo
$(mc68)_ciface.lo: deps90.lo

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

CLEANFILES = $(HSL_SHARED_LIB) $(DEFFILE) *.mod

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
