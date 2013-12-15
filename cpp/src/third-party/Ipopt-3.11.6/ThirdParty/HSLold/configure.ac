# Copyright (C) 2006, 2009 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: configure.ac 2927 2013-09-12 12:55:38Z stefan $

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([ThirdPartyHSL],[1.3.6],[http://projects.coin-or.org/BuildTools/newticket])

AC_COPYRIGHT([
Copyright 2006, 2009 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
# Here we don't know yet which files are actually provided...
AC_CONFIG_SRCDIR(configure.ac)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

AC_COIN_PROJECTDIR_INIT(COINHSL,4:6:3)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# turn off warnings and pedantic mode since they could make trouble on some platforms,
# and we couldn't fix them anyway
coin_skip_warn_fflags=yes

# provide flag to allow compilation of a stand-alone shared library
AC_ARG_ENABLE([loadable-library],
[AC_HELP_STRING([--enable-loadable-library],
                [compile stand-alone shared library for delayed loading])])

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(HSL)

# Get the name of the C compiler (for a dummy program)
AC_COIN_PROG_CC

# Get the name of the Fortran compiler and appropriate compiler options
AC_COIN_PROG_F77

# Check if pkg-config is available (to setup addlibs file)
AC_COIN_HAS_PKGCONFIG

# Get the linker flags (FLIBS) that are required for linking against this blas library
AC_COIN_F77_WRAPPERS

# Add FLIBS to HSL_PCLIBS, so that they get into the .pc files
HSL_PCLIBS="$HSL_PCLIBS $FLIBS"

# Need egrep later on
AC_PROG_EGREP

# Also set up the FC compiler variable, required to compile .f90
# files.  We use the same compile and flags as F77 and FFLAGS, unless
# the user has set FC and FCFLAGS explicitly.
if test "$FC" = ""; then
  FC="$F77"
fi
if test "$FCFLAGS" = ""; then
  FCFLAGS="$FFLAGS"
fi
AC_PROG_FC

# Initialize autotools

if test x$enable_loadable_library = xyes; then
  flag=force_shared
fi
AC_COIN_INIT_AUTO_TOOLS([$flag])

#############################################################################
#                    Find out which files are available                     #
#############################################################################

AC_MSG_CHECKING([whether we can compile Fortran 90 code])
AC_LANG_PUSH(Fortran 77)
AC_COMPILE_IFELSE([        MODULE BLA
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BLUBB
      END MODULE BLA],[have_f90=yes ; rm -f bla.mod],[have_f90=no])
AC_MSG_RESULT([$have_f90])
AC_LANG_POP(Fortran 77)

coin_has_something=no

AC_DEFUN([MY_CHECK_HSLFILE],
[AC_COIN_CHECK_FILE([$srcdir/$1.f],
                    [hsl_files="$hsl_files $1.f"
                     AC_DEFINE(COINHSL_HAS_$2, [1], [Define to 1 if $2 is available.])
                     coin_has_something=yes
                     coin_has_$1=yes],
                    [coin_has_$1=no])
AM_CONDITIONAL([COIN_HAS_$2],[test $coin_has_$1 = yes])])

AC_DEFUN([MY_CHECK_HSLFILE_F90],
[if test $have_f90 = yes; then
   AC_COIN_CHECK_FILE([$srcdir/$1.f90],
                      [hsl_files="$hsl_files $1.f"
                       AC_DEFINE(COINHSL_HAS_$2, [1], [Define to 1 if $2 is available.])
                       coin_has_something=yes
                       coin_has_$1=yes],
                      [coin_has_$1=no])
else
  coin_has_$1=no
fi
AM_CONDITIONAL([COIN_HAS_$2],[test $coin_has_$1 = yes])])

MY_CHECK_HSLFILE(ma27ad,MA27)
if test $coin_has_ma27ad = yes; then
  AC_MSG_CHECKING([whether MA27 is threadsafe])
  if $EGREP MA27UD $srcdir/ma27ad.f >/dev/null 2>&1; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Please obtain the threadsafe (newer) version of MA27])
  fi
fi

AC_COIN_CHECK_PACKAGE(Metis, [coinmetis], HSL)
if test $coin_has_metis = yes ; then
  AC_DEFINE(COINHSL_HAS_METIS, [1], Define to 1 if METIS is available)
fi

MY_CHECK_HSLFILE(ma57ad,MA57)
if test $coin_has_ma57ad = yes; then
  if test $coin_has_metis != yes ; then
    AC_MSG_ERROR([MA57 requires METIS, which is not found.])
  fi
  AC_COIN_CHECK_PACKAGE_BLAS(HSL)
  if test $coin_has_blas != yes ; then
    AC_MSG_ERROR([MA57 requires BLAS, which is not found.])
  fi
else
  AM_CONDITIONAL(COIN_HAS_METIS, [test 0 = 1])
  AM_CONDITIONAL(COIN_HAS_BLAS, [test 0 = 1])
fi

MY_CHECK_HSLFILE(mc19ad,MC19)
MY_CHECK_HSLFILE(ma28ad,MA28)
MY_CHECK_HSLFILE_F90(hsl_ma77d,MA77)

#############################################################################
#                        Setup build of loadable library                    #
#############################################################################

HSL_SHARED_LIB=
if test x$enable_loadable_library = xyes; then
  case $build in
    *-cygwin* | *-mingw*)
      HSL_SHARED_LIB=libhsl.dll
      F77_SHARED_CMD='$(F77) -shared $(FFLAGS)'
      case $F77 in
        compile_f2c*)
          EXPORT_SYMBOLS=
          if test $coin_has_ma27ad = yes; then
	    EXPORT_SYMBOLS="$EXPORT_SYMBOLS ma27id_ ma27ad_ ma27bd_ ma27cd_"
          fi
          if test $coin_has_ma28ad = yes; then
	    EXPORT_SYMBOLS="$EXPORT_SYMBOLS ma28ad_"
          fi
          if test $coin_has_ma57ad = yes; then
	    EXPORT_SYMBOLS="$EXPORT_SYMBOLS ma57id_ ma57ad_ ma57bd_ ma57cd_ ma57ed_"
          fi
          if test $coin_has_mc19ad = yes; then
	    EXPORT_SYMBOLS="$EXPORT_SYMBOLS mc19ad_"
          fi
          F77_SHARED_CMD='$(F77) -shared $(FFLAGS) -def:$(DEFFILE)'
	  DEFFILE='Export.def'
          ;;
        ifort* | IFORT*)
	  EXPORT_SYMBOLS=
          if test $coin_has_ma27ad = yes; then
	    EXPORT_SYMBOLS="$EXPORT_SYMBOLS MA27ID MA27AD MA27BD MA27CD"
          fi
          if test $coin_has_ma28ad = yes; then
	    EXPORT_SYMBOLS="$EXPORT_SYMBOLS MA28AD"
          fi
          if test $coin_has_ma57ad = yes; then
	    EXPORT_SYMBOLS="$EXPORT_SYMBOLS MA57ID MA57AD MA57BD MA57CD MA57ED"
          fi
          if test $coin_has_mc19ad = yes; then
	    EXPORT_SYMBOLS="$EXPORT_SYMBOLS MC19AD"
          fi
	  DEFFILE='Export.def'
          F77_SHARED_CMD='$(F77) $(FFLAGS)'
	  F77_DLL_LINKFLAGS='-link -dll -def:$(DEFFILE)'
          ;;
        *)
          F77_SHARED_CMD='$(F77) -shared $(FFLAGS)'
          ;;
      esac
      ;;
  esac
fi
AC_SUBST(F77_SHARED_CMD)
AC_SUBST(F77_DLL_LINKFLAGS)
AC_SUBST(HSL_SHARED_LIB)
AC_SUBST(EXPORT_SYMBOLS)
AC_SUBST(DEFFILE)
AM_CONDITIONAL([INSTALL_LOADABLE_LIBRARY],[test x$flag = xforce_shared])
AM_CONDITIONAL([NEED_DEFFILE], [test ! x"$DEFFILE" = x])

AM_CONDITIONAL([EXTRA_SHARED_COMMAND],[test ! x"$F77_SHARED_CMD" = x])
#AM_CONDITIONAL([EXTRA_SHARED_COMMAND],[test "$F77_SHARED_CMD" = x])

#############################################################################
#                              Write the output                             #
#############################################################################

AC_CONFIG_FILES([Makefile])

# we build the library and setup the .pc files only if there is at least one source code
# if there is no source code, then .pc files should not be setup, so for others this is as if we do not exist
if test $coin_has_something = yes ; then
  AC_CONFIG_FILES([coinhsl.pc coinhsl-uninstalled.pc])
  AC_CONFIG_HEADER(config.h config_coinhsl.h)
else
  # in case something is left over from a previous build
  rm -f coinhsl.pc coinhsl-uninstalled.pc
fi
AM_CONDITIONAL(BUILD_COINHSL, [test $coin_has_something = yes])

AC_COIN_FINALIZE
