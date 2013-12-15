# Copyright (C) 2006, 2009 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: configure.ac 3092 2013-11-12 15:40:26Z stefan $

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([ThirdPartyHSL],[1.4.4],[http://projects.coin-or.org/BuildTools/newticket])

AC_COPYRIGHT([
Copyright 2006, 2009 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package BuildTools which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR([coinhsl/common/deps.f])

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

AC_COIN_PROJECTDIR_INIT(COINHSL,5:4:4)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

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

# recent autotools versions can recognize OpenMP with the following macro
# this would be useful when using MA86 and MA97 
# AC_OPENMP

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
      PUBLIC TEST_TYPE
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BLUBB
      TYPE TEST_TYPE
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TTVAR
      END TYPE TEST_TYPE
      END MODULE BLA],[have_f90=yes ; rm -f bla.mod],[have_f90=no])
AC_MSG_RESULT([$have_f90])
AC_LANG_POP(Fortran 77)

coin_has_something=no
hsl2013=no

AC_DEFUN([MY_CHECK_HSLFILE],
[AC_MSG_CHECKING([for m4_toupper(m4_bpatsubst($1, hsl_))])
AC_COIN_CHECK_FILE([$srcdir/coinhsl/$1/$1d.f],
                    [AC_DEFINE(m4_toupper(COINHSL_HAS_$1), [1], Define to 1 if m4_toupper(m4_bpatsubst($1, hsl_)) is available.)
                     coin_has_something=yes
                     coin_has_$1=yes],
                    [coin_has_$1=no])
AC_MSG_RESULT($coin_has_$1)
AM_CONDITIONAL(m4_toupper(COIN_HAS_$1),[test $coin_has_$1 = yes])
])

AC_DEFUN([MY_CHECK_HSLFILE_F90],
[if test $have_f90 = yes; then
   AC_MSG_CHECKING([for m4_toupper(m4_bpatsubst($1, hsl_))])
   AC_COIN_CHECK_FILE([$srcdir/coinhsl/$1/$1d.f90],
                      [AC_DEFINE([COINHSL_HAS_]m4_toupper(m4_bpatsubst($1, hsl_)), [1], Define to 1 if m4_toupper(m4_bpatsubst($1, hsl_)) is available.)
                       coin_has_something=yes
                       coin_has_$1=yes],
                      [coin_has_$1=no])
   AC_MSG_RESULT($coin_has_$1)
else
  coin_has_$1=no
fi
AM_CONDITIONAL(m4_toupper(COIN_HAS_$1),[test $coin_has_$1 = yes])
])

MY_CHECK_HSLFILE(mc19)
MY_CHECK_HSLFILE(ma27)
MY_CHECK_HSLFILE(ma28)
MY_CHECK_HSLFILE(ma57)
MY_CHECK_HSLFILE_F90(hsl_ma77)
MY_CHECK_HSLFILE_F90(hsl_ma86)
MY_CHECK_HSLFILE_F90(hsl_ma97)

if test $have_f90 = yes; then
  AC_MSG_CHECKING([for C interface to MC68])
  AC_COIN_CHECK_FILE([$srcdir/coinhsl/hsl_mc68/C/hsl_mc68i_ciface.f90],
                     [AC_DEFINE(COINHSL_HAS_MC68, [1], Define to 1 if MC68 is available.)
                      coin_has_something=yes
                      grep -q mc68_default_control_i $srcdir/coinhsl/hsl_mc68/C/hsl_mc68i_ciface.f90 && hsl2013=yes
                      coin_has_hsl_mc68=yes],
                     [coin_has_hsl_mc68=no])
  AC_MSG_RESULT($coin_has_hsl_mc68)
else
  coin_has_hsl_mc68=no
fi
AM_CONDITIONAL(COIN_HAS_HSL_MC68,[test $coin_has_hsl_mc68 = yes])

# the F90 codes depend on deps90.f90
if test $coin_has_hsl_mc68$coin_has_hsl_ma77$coin_has_hsl_ma86$coin_has_hsl_ma97 != nononono ; then
  AC_MSG_CHECKING([for deps90.f90])
  AC_COIN_CHECK_FILE([$srcdir/coinhsl/common/deps90.f90],
                     [coin_has_depsf90=yes],
                     [coin_has_depsf90=no])
  AC_MSG_RESULT($coin_has_depsf90)
  if test $coin_has_depsf90 = no; then
    AC_MSG_ERROR([Required file common/deps90.f90 not found])
  fi
else
  coin_has_depsf90=no
fi
AM_CONDITIONAL(COIN_HAS_DEPSF90,[test $coin_has_depsf90 = yes])

# MA57, MA77, and MA86 require BLAS
if test $coin_has_ma57$coin_has_hsl_ma77$coin_has_hsl_ma86 != nonono ; then
  AC_COIN_CHECK_PACKAGE_BLAS(HSL)
  if test $coin_has_blas != yes ; then
    AC_MSG_ERROR([Required package BLAS not found.])
  fi
else
  AM_CONDITIONAL(COIN_HAS_BLAS, [test 0 = 1])
fi

# MA57, MA77, MA97, and MC68 can use Metis
# if Metis is not available, then we compile a dummy metis.f that defines the Metis function symbol
if test $coin_has_ma57$coin_has_hsl_ma77$coin_has_hsl_ma97$coin_has_hsl_mc68 != nononono ; then
  AC_COIN_CHECK_PACKAGE(Metis, [coinmetis], HSL)
  if test $coin_has_metis = yes ; then
    AC_DEFINE(COINHSL_HAS_METIS, [1], Define to 1 if METIS is available)
  fi
else
  # do as if we had metis, so metis.f is not compiled
  AM_CONDITIONAL(COIN_HAS_METIS, [test 1 = 1])
fi

# MC68 and MA77 changed the interface from HSL 2012 to HSL 2013,
# so we check whether we have the 2013 version (for mc68, we did this test above)
# and set a compiletime define
if test $coin_has_hsl_ma77 = yes ; then
  grep -q ma77_solve_fredholm $srcdir/coinhsl/hsl_ma77/hsl_ma77d.f90 && hsl2013=yes
fi
if test $hsl2013 = yes ; then
  AC_DEFINE(COINHSL_HSL2013, [1], Define to 1 if HSL library is from 2013)
  AC_MSG_NOTICE([HSL source recognized as HSL 2013])
fi

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
          if test $coin_has_ma27 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS ma27ad_ ma27bd_ ma27cd_ ma27id_"
          fi
          if test $coin_has_ma28 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS ma28ad_"
          fi
          if test $coin_has_ma57 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS ma57ad_ ma57bd_ ma57cd_ ma57ed_ ma57id_"
          fi
          if test $coin_has_mc19 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS mc19ad_"
          fi
          if test $coin_has_hsl_ma77 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS ma77_default_control_ ma77_open_nelt_ ma77_open_ ma77_input_vars_ ma77_input_reals_ ma77_analyse_ ma77_factor_ ma77_factor_solve_ ma77_solve_ ma77_resid_ ma77_scale_ ma77_enquire_posdef_ ma77_enquire_indef_ ma77_alter_ ma77_restart_ ma77_finalise_"            
          fi
          if test $coin_has_hsl_ma86 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS ma86_default_control_ ma86_analyse_ ma86_factor_ ma86_factor_solve_ ma86_solve_ ma86_finalise_"
          fi
          if test $coin_has_hsl_ma97 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS ma97_default_control_ ma97_analyse_ ma97_factor_ ma97_factor_solve_ ma97_solve_ ma97_finalise_ ma97_free_akeep_"
          fi
          if test $coin_has_hsl_mc68 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS mc68_default_control_ mc68_order_"
          fi
          F77_SHARED_CMD='$(F77) -shared $(FFLAGS) -def:$(DEFFILE)'
          DEFFILE='Export.def'
          ;;
        ifort* | IFORT*)
          EXPORT_SYMBOLS=
          if test $coin_has_ma27 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS MA27AD MA27BD MA27CD MA27ID"
          fi
          if test $coin_has_ma28 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS MA28AD"
          fi
          if test $coin_has_ma57 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS MA57AD MA57BD MA57CD MA57ED MA57ID"
          fi
          if test $coin_has_hsl_ma77 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS MA77_DEFAULT_CONTROL MA77_OPEN_NELT MA77_OPEN MA77_INPUT_VARS MA77_INPUT_REALS MA77_ANALYSE MA77_FACTOR MA77_FACTOR_SOLVE MA77_SOLVE MA77_RESID MA77_SCALE MA77_ENQUIRE_POSDEF MA77_ENQUIRE_INDEF MA77_ALTER MA77_RESTART MA77_FINALISE"            
          fi
          if test $coin_has_mc19 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS MC19AD"
          fi
          if test $coin_has_hsl_ma86 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS MA86_DEFAULT_CONTROL MA86_ANALAYSE MA86_FACTOR MA86_FACTOR_SOLVE MA86_SOLVE MA86_FINALISE"
          fi
          if test $coin_has_hsl_ma97 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS MA97_DEFAULT_CONTROL MA97_ANALAYSE MA97_FACTOR MA97_FACTOR_SOLVE MA97_SOLVE MA97_FINALISE MA97_FREE_AKEEP"
          fi
          if test $coin_has_hsl_mc68 = yes; then
            EXPORT_SYMBOLS="$EXPORT_SYMBOLS MC68_DEFAULT_CONTROL MC68_ORDER"
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
