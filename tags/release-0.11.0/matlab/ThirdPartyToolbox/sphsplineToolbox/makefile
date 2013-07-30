#		$Id: makefile,v 1.4 2012/03/14 23:31:10 myself Exp $
#
#	Makefile for sphspline: spherical surface spline mex functions
#	build in C.  Suitable for Matlab or Octave; the choice is made
#	during the GMT installation (rerun GMT's configure script if
#	you need to switch as only one mex interface can be built.)
#	Use with GMT 4.5.2 or later
#
#	Paul Wessel, May 2008.  Revised 2010.
#	Ref: Wessel, P. and J. M .Becker, 2008, "Interpolation using a
#	  generalized Green's function for a spherical surface spline
#	  in tension", Geophys. J. Int., doi:10.1111/j.1365-246X.2008.03829.x
#
# Note that if you have write permission in the Matlab/Octave install dir
# then that is where install will place the files.  Otherwise, it
# will put them in $(bindir).
#----------------------------------------------------------------------------
#	The rest should should be ok. NO EDITING BEYOND HERE, unless your
#	Matlab setup was not auto-detected by the configure file
#	If so, please let us know what you needed to change.
#----------------------------------------------------------------------------

GMTSRCDIR	= $(GMTHOME)/src/
include $(GMTSRCDIR)config.mk
include $(GMTSRCDIR)common.mk

INCLUDES	= -I$(srcdir)  $(NETCDF_INC)

.SUFFIXES:	.m .$(MEX_EXT)

ALLLIB     	= $(GMT_LIB) $(NETCDF_LIB) $(GDAL_LIB) $(MEX_LIB)
FLAGS		= -I$(srcdir) $(MEX_INC) $(NETCDF_INC)

PROGS_M		= SSSTc.m SSSTgradc.m 
PROGS_MEX	= SSSTc.$(MEX_EXT) SSSTgradc.$(MEX_EXT) WB_lookup.$(MEX_EXT)

#-------------------------------------------------------------------------------
#	software targets
#-------------------------------------------------------------------------------

all:		$(PROGS_MEX)

install:	all
		$(INSTALL) -d $(MEX_XDIR)
		if [ -w $(MEX_XDIR) ] ; then \
			$(INSTALL) $(PROGS_MEX) $(MEX_XDIR); \
		else \
			$(INSTALL) $(PROGS_MEX) $(bindir); \
		fi
		$(INSTALL) -d $(MEX_MDIR)
		if [ -w $(MEX_MDIR) ] ; then \
			$(INSTALL) $(PROGS_M) $(MEX_MDIR); \
		else \
			$(INSTALL) $(PROGS_M) $(bindir); \
		fi

uninstall:
		if [ -w $(MEX_XDIR) ] ; then \
			cd $(MEX_XDIR); \
		else \
			cd $(bindir); \
		fi; \rm -f $(PROGS_MEX)
		if [ -w $(MEX_MDIR) ] ; then \
			cd $(MEX_MDIR); \
		else \
			cd $(bindir); \
		fi; \rm -f $(PROGS_M)

install-man uninstall-man manpages:	;

spotless::	clean

clean:
		rm -f .gmt* *.o *% .mexrc.sh $(PROGS_MEX)

#-------------------------------------------------------------------------------
#	program rules
#-------------------------------------------------------------------------------

SSSTc.$(MEX_EXT):	$(GMT_H) $(LIBGMT) SSSTc.c SSSTc.m
		$(MEX_BLD) $(FLAGS) SSSTc.c $(ALLLIB)

SSSTgradc.$(MEX_EXT):	$(GMT_H) $(LIBGMT) SSSTgradc.c SSSTgradc.m
		$(MEX_BLD) $(FLAGS) SSSTgradc.c $(ALLLIB) 

WB_lookup.$(MEX_EXT):	$(GMT_H) $(LIBGMT) WB_lookup.c WB_lookup.m
		$(MEX_BLD) $(FLAGS) WB_lookup.c $(ALLLIB)
