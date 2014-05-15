#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the class library                   *#
#*       SoPlex --- the Sequential object-oriented simPlex.                  *#
#*                                                                           *#
#*    Copyright (C) 1996-2012 Konrad-Zuse-Zentrum                            *#
#*                            fuer Informationstechnik Berlin                *#
#*                                                                           *#
#*  SoPlex is distributed under the terms of the ZIB Academic Licence.       *#
#*                                                                           *#
#*  You should have received a copy of the ZIB Academic License              *#
#*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  *#
#*                                                                           *#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#

#@file    Makefile
#@brief   SoPlex Makefile
#@author  Thorsten Koch

#-----------------------------------------------------------------------------
# detect host architecture
#-----------------------------------------------------------------------------
include make/make.detecthost


#-----------------------------------------------------------------------------
VERSION		:=	1.7.2

VERBOSE		=	false
SHARED		=	false
OPT		=	opt
STATICLIBEXT	=	a
SHAREDLIBEXT	=	so
LIBEXT		=	$(STATICLIBEXT)
EXEEXTENSION	=	
TEST		=	quick
ALGO		=	1 2 3 4
LIMIT		=	#

#these variables are needed for cluster runs
TIME		=	3600
MEM		=	6144
CONTINUE	=	false

INSTALLDIR	=	#

#will this be compiled for PARASCIP? (disables output because it uses global variables)
PARASCIP	=	false

GMP		=	false
ZLIB		=	true

COMP		=	gnu
CXX		=	g++
CXX_c		=	-c # the trailing space is important
CXX_o		=	-o # the trailing space is important
LINKCXX		=	$(CXX)
LINKCXX_L	=	-L
LINKCXX_l	=	-l
LINKCXX_o	=	-o # the trailing space is important
LINKLIBSUFFIX	=
DCXX		=	$(CXX)
LINT		=	flexelint
AR		=	ar
AR_o		=
RANLIB		=	ranlib
DOXY		=	doxygen
VALGRIND	=	valgrind
LN_s		=	ln -s

LIBBUILD	=	$(AR)
LIBBUILD_o	=	$(AR_o)
LIBBUILDFLAGS	=       $(ARFLAGS)

CPPFLAGS	=	-Isrc
CXXFLAGS	=	
BINOFLAGS	=	
LIBOFLAGS	=	
LDFLAGS		=	
ARFLAGS		=	cr
DFLAGS		=	-MM
VFLAGS		=	--tool=memcheck --leak-check=yes --show-reachable=yes #--gen-suppressions=yes

GMP_FLAGS	=
GMP_LDFLAGS	=	-lgmpxx -lgmp

SOPLEXDIR	=	$(realpath .)
SRCDIR		=	src
BINDIR		=	bin
LIBDIR		=	lib
INCLUDEDIR	=	include
NAME		=	soplex
LIBOBJ		= 	changesoplex.o clufactor.o didxset.o \
			dsvector.o dvector.o dvector_exact.o enter.o \
			idxset.o leave.o lpcolset.o lprowset.o \
			lprow.o mpqreal.o mpsinput.o nameset.o \
			slufactor.o soplex.o \
			spxbasis.o spxbounds.o spxboundflippingrt.o spxchangebasis.o \
			spxequilisc.o spxdantzigpr.o spxdefaultrt.o \
			spxdefines.o spxdesc.o spxdevexpr.o \
			spxfastrt.o spxfileio.o spxgeometsc.o spxgithash.o\
			spxharrisrt.o spxhybridpr.o spxid.o spxio.o \
			spxlp.o spxlpfread.o spxmainsm.o spxmpsread.o \
			spxmpswrite.o spxlpfwrite.o \
			spxout.o spxparmultpr.o spxquality.o \
			spxscaler.o spxshift.o spxsolver.o spxsolve.o \
			spxstarter.o spxsteeppr.o spxsumst.o spxvecs.o \
			spxvectorst.o spxweightpr.o spxweightst.o spxwritestate.o \
			ssvector.o svector.o svset.o timer.o \
			tracemethod.o unitvector.o updatevector.o \
			vector.o vector_exact.o \
			gzstream.o
BINOBJ		=	soplexmain.o
EXAMPLEOBJ	=	simpleexample.o
REPOSIT		=	# template repository, explicitly empty  #spxproof.o 

BASE		=	$(OSTYPE).$(ARCH).$(COMP).$(OPT)

LASTSETTINGS	=	$(OBJDIR)/make.lastsettings


#------------------------------------------------------------------------------
#--- NOTHING TO CHANGE FROM HERE ON -------------------------------------------
#------------------------------------------------------------------------------

GCCWARN		=	-Wall -W -Wpointer-arith -Wno-unknown-pragmas \
			-Wcast-align -Wwrite-strings -Wconversion \
			-Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder \
			-Woverloaded-virtual -Wsign-promo -Wsynth -Wundef \
			-Wcast-qual -Wold-style-cast -Wshadow 
#			-Weffc++ -Wredundant-decls    
# gcc 2.xx -Wmissing-declarations -Wbad-function-cast 

#GCCWARN =
#-----------------------------------------------------------------------------
include make/make.$(BASE)
-include make/local/make.$(HOSTNAME)
-include make/local/make.$(HOSTNAME).$(COMP)
-include make/local/make.$(HOSTNAME).$(COMP).$(OPT)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# SHARED Libaries
#-----------------------------------------------------------------------------

ifeq ($(SHARED),true)
CPPFLAGS	+=	-fPIC
LIBEXT		=	$(SHAREDLIBEXT)
LIBBUILD	=	$(LINKCXX)
LIBBUILDFLAGS	+=      -shared
LIBBUILD_o	= 	-o # the trailing space is important
ARFLAGS		=
RANLIB		=
endif

CXXFLAGS	+=	$(USRCXXFLAGS)
LDFLAGS		+=	$(USRLDFLAGS)
ARFLAGS		+=	$(USRARFLAGS)
DFLAGS		+=	$(USRDFLAGS)

#-----------------------------------------------------------------------------
# PARASCIP
#-----------------------------------------------------------------------------

ifeq ($(PARASCIP),true)
CPPFLAGS	+=	-DDISABLE_VERBOSITY
endif

#-----------------------------------------------------------------------------

BINNAME		=	$(NAME)-$(VERSION).$(BASE)
EXAMPLENAME	=	simpleexample.$(BASE)
LIBNAME		=	$(NAME)-$(VERSION).$(BASE)
BINFILE		=	$(BINDIR)/$(BINNAME)$(EXEEXTENSION)
EXAMPLEFILE	=	$(BINDIR)/$(EXAMPLENAME)
LIBFILE		=	$(LIBDIR)/lib$(LIBNAME).$(LIBEXT)
LIBSHORTLINK	=	$(LIBDIR)/lib$(NAME).$(LIBEXT)
LIBLINK		=	$(LIBDIR)/lib$(NAME).$(BASE).$(LIBEXT)
BINLINK		=	$(BINDIR)/$(NAME).$(BASE)$(EXEEXTENSION)
BINSHORTLINK	=	$(BINDIR)/$(NAME)$(EXEEXTENSION)
DEPEND		=	src/depend

# potential valgrind suppression file name
VSUPPNAME	= 	$(OSTYPE).$(ARCH).$(COMP).supp

OBJDIR		=	obj/O.$(BASE)
BINOBJDIR	=	$(OBJDIR)/bin
LIBOBJDIR	=	$(OBJDIR)/lib
BINOBJFILES	=	$(addprefix $(BINOBJDIR)/,$(BINOBJ))
EXAMPLEOBJFILES	=	$(addprefix $(BINOBJDIR)/,$(EXAMPLEOBJ))
LIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(LIBOBJ))
BINSRC		=	$(addprefix $(SRCDIR)/,$(BINOBJ:.o=.cpp))
EXAMPLESRC	=	$(addprefix $(SRCDIR)/,$(EXAMPLEOBJ:.o=.cpp))
LIBSRC		=	$(addprefix $(SRCDIR)/,$(LIBOBJ:.o=.cpp))
LIBSRCHEADER	=	$(addprefix $(SRCDIR)/,$(LIBOBJ:.o=.h))

GMPDEP		:=	$(SRCDIR)/depend.gmp
GMPSRC		:=	$(shell cat $(GMPDEP))
ifeq ($(GMP),true)
CPPFLAGS	+=	-DSOPLEX_WITH_GMP $(GMP_FLAGS)
LDFLAGS		+=	$(GMP_LDFLAGS)
endif

ZLIBDEP		:=	$(SRCDIR)/depend.zlib
ZLIBSRC		:=	$(shell cat $(ZLIBDEP))
ifeq ($(ZLIB_LDFLAGS),)
ZLIB		=	false
endif
ifeq ($(ZLIB),true)
CPPFLAGS	+=	-DWITH_ZLIB $(ZLIB_FLAGS)
LDFLAGS		+=	$(ZLIB_LDFLAGS)
endif


ifeq ($(VERBOSE),false)
.SILENT:	$(LIBLINK) $(LIBSHORTLINK) $(BINLINK) $(BINSHORTLINK) $(BINFILE) $(EXAMPLEFILE) $(EXAMPLEOBJFILES) $(LIBFILE) $(BINOBJFILES) $(LIBOBJFILES)
endif

all:		githash $(LIBFILE) $(BINFILE) $(LIBLINK) $(LIBSHORTLINK) $(BINLINK) $(BINSHORTLINK)

simpleexample:	$(LIBFILE) $(EXAMPLEFILE) $(LIBLINK) $(LIBSHORTLINK)

# include install targets
-include make/make.install

$(LIBLINK) $(LIBSHORTLINK):	$(LIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(LIBFILE)) $(notdir $@)

$(BINLINK) $(BINSHORTLINK):	$(BINFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(BINFILE)) $(notdir $@)

$(BINFILE):	$(BINDIR) $(BINOBJDIR) $(LIBOBJFILES) $(BINOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(BINOBJFILES) $(LIBOBJFILES) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(EXAMPLEFILE):	$(BINDIR) $(EXAMPLEOBJDIR) $(LIBOBJFILES) $(EXAMPLEOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(EXAMPLEOBJFILES) $(LIBOBJFILES) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(LIBFILE):	$(LIBDIR) $(LIBOBJDIR) touchexternal $(LIBOBJFILES)
		@echo "-> generating library $@"
		-rm -f $(LIBFILE)
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(LIBOBJFILES) $(REPOSIT)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

-include make/local/make.targets

.PHONY: lint
lint:		$(BINSRC) $(LIBSRC)
		$(LINT) lint/$(NAME).lnt -os\(lint.out\) \
		$(CPPFLAGS) -UNDEBUG $^

.PHONY: doc
doc:		
		cd doc; $(DOXY) $(NAME).dxy

.PHONY: test
test:		check

.PHONY: check
check:		#$(BINFILE)
		cd check; ./check.sh $(TEST).test ../$(BINFILE) '$(ALGO)' $(LIMIT)

valgrind-check:	$(BINFILE)
		cd check; \
		./valgrind.sh $(TEST).test ../$(BINFILE) '$(ALGO)' '$(LIMIT)' \
		"$(VALGRIND) $(VFLAGS)" $(VSUPPNAME)

memory_exception_test: $(BINFILE)
		cd check; \
		./exception.sh $(TEST).test ../$(BINFILE) '$(ALGO)' '$(LIMIT)' \
		"$(VALGRIND) $(VFLAGS)" $(VSUPPNAME)

.PHONY: cleanbin
cleanbin:       $(BINDIR)
		@echo "remove binary $(BINFILE)"
		@-rm -f $(BINFILE) $(BINLINK) $(BINSHORTLINK)

.PHONY: cleanlib
cleanlib:       $(LIBDIR)
		@echo "remove library $(LIBFILE)" 
		@-rm -f $(LIBFILE) $(LIBLINK) $(LIBSHORTLINK)

.PHONY: clean
clean:          cleanlib cleanbin $(LIBOBJDIR) $(BINOBJDIR) $(OBJDIR)
		@echo "remove objective files" 
ifneq ($(LIBOBJDIR),)
		@-rm -f $(LIBOBJDIR)/*.o && rmdir $(LIBOBJDIR)
endif
ifneq ($(BINOBJDIR),)
		@-rm -f $(BINOBJDIR)/*.o && rmdir $(BINOBJDIR)
endif
ifneq ($(OBJDIR),)
		@-rm -f $(LASTSETTINGS)
		@-rmdir $(OBJDIR)
endif
		@-rm -f $(EXAMPLEFILE)

vimtags:
		-ctags -o TAGS src/*.cpp src/*.h

etags:
		-ctags -e -o TAGS src/*.cpp src/*.h

$(OBJDIR):	
		@-mkdir -p $(OBJDIR)

$(BINOBJDIR):	$(OBJDIR)
		@-mkdir -p $(BINOBJDIR)

$(LIBOBJDIR):	$(OBJDIR)
		@-mkdir -p $(LIBOBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

$(LIBDIR):
		@-mkdir -p $(LIBDIR)

.PHONY: depend
depend:
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) \
		$(BINSRC:.o=.cpp) \
		| sed '\''s|^\([0-9A-Za-z_]\{1,\}\)\.o|$$\(BINOBJDIR\)/\1.o|g'\'' \
		>$(DEPEND)'
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) \
		$(EXAMPLESRC:.o=.cpp) \
		| sed '\''s|^\([0-9A-Za-z_]\{1,\}\)\.o|$$\(BINOBJDIR\)/\1.o|g'\'' \
		>>$(DEPEND)'
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) \
		$(LIBSRC:.o=.cpp) \
		| sed '\''s|^\([0-9A-Za-z_]\{1,\}\)\.o|$$\(LIBOBJDIR\)/\1.o|g'\'' \
		>>$(DEPEND)'
		@echo `grep -l "SOPLEX_WITH_GMP" $(SRCDIR)/*` >$(GMPDEP)
		@echo `grep -l "WITH_ZLIB" $(SRCDIR)/*` >$(ZLIBDEP)

-include	$(DEPEND)

$(BINOBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@-mkdir -p $(BINOBJDIR)
		@echo "-> compiling $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BINOFLAGS) $(CXX_c)$< $(CXX_o)$@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@-mkdir -p $(LIBOBJDIR)
		@echo "-> compiling $@"
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LIBOFLAGS) $(CXX_c)$< $(CXX_o)$@


-include $(LASTSETTINGS)

.PHONY: touchexternal
touchexternal:	$(GMPDEP) $(ZLIBDEP)
ifneq ($(GMP),$(LAST_GMP))
		@-touch $(GMPSRC)
endif
ifneq ($(ZLIB),$(LAST_ZLIB))
		@-touch $(ZLIBSRC)
endif
ifneq ($(SHARED),$(LAST_SHARED))
		@-touch $(LIBSRC)
		@-touch $(BINSRC)
endif
ifneq ($(USRCXXFLAGS),$(LAST_USRCXXFLAGS))
		@-touch $(LIBSRC)
		@-touch $(BINSRC)
endif
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_GMP=$(GMP)" >> $(LASTSETTINGS)
		@echo "LAST_ZLIB=$(ZLIB)" >> $(LASTSETTINGS)
		@echo "LAST_SHARED=$(SHARED)" >> $(LASTSETTINGS)
		@echo "LAST_USRCXXFLAGS=$(USRCXXFLAGS)" >> $(LASTSETTINGS)

-include make/local/make.detectgithash
# this empty target is needed for the SoPlex release versions
githash::	# do not remove the double-colon

# --- EOF ---------------------------------------------------------------------
