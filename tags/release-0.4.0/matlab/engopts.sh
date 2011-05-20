#
# engopts.sh	Shell script for configuring engine standalone applications.
#               These options were tested with the specified compiler.
#
# usage:        Do not call this file directly; it is sourced by the
#               mbuild shell script.  Modify only if you don't like the
#               defaults after running mbuild.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# Note: For the version of system compiler supported with this release,
#       refer to the Supported and Compatible Compiler List at:
#       http://www.mathworks.com/support/compilers/current_release/
#
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_SA_OPT: Template Options file for building standalone engine applications
#
# Copyright 1984-2008 The MathWorks, Inc.
# $Revision: 1.30.4.13 $  $Date: 2008/11/04 19:40:05 $
#----------------------------------------------------------------------------
#
    if [ "$TMW_ROOT" = "" ]; then
	TMW_ROOT="$MATLAB"
    fi
    MFLAGS="-I$TMW_ROOT/extern/include"
    MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat"
    MCXXFLAGS="-I$TMW_ROOT/extern/include/cpp $MFLAGS"
    MBAFLAGS="-I`pwd`/../cpp/src/third-party/mba/include/"
    MBALIBS="-L`pwd`/PointsToolbox -lMBA"
    MCXXLIBS="$MBALIBS $MLIBS"
    LDEXTENSION=''
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The mex script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        glnx86)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath=$TMW_ROOT/bin/$Arch -Wl,-rpath=`pwd`/PointsToolbox"
            CC='gcc-4.3'
            CFLAGS='-ansi -D_GNU_SOURCE'
            CFLAGS="$CFLAGS -fexceptions"
            CFLAGS="$CFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CFLAGS="$CFLAGS $MFLAGS"
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
            CLIBS="$CLIBS -lstdc++"
#           
            CXX='g++-4.3'
            CXXFLAGS='-ansi -D_GNU_SOURCE'
            CXXFLAGS="$CXXFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CXXFLAGS="$CXXFLAGS $MCXXFLAGS $MBAFLAGS -DGLNX86 -DGCC"
            CXXLIBS="$RPATH $MCXXLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
#
            FC='g95'
            FFLAGS='-fexceptions'
            FFLAGS="$FFLAGS $MFLAGS"
            FLIBS="$RPATH $MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-pthread -shared -Wl,--version-script,$TMW_ROOT/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined"
	    LDFLAGS="$LDFLAGS"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
	    LDEXTENSION='.mexglx'
#
            POSTLINK_CMDS=""
#----------------------------------------------------------------------------
            ;;
        glnxa64)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath=$TMW_ROOT/bin/$Arch -Wl,-rpath=`pwd`/PointsToolbox"
            CC='gcc-4.3'
            CFLAGS='-ansi -D_GNU_SOURCE'
            CFLAGS="$CFLAGS -fexceptions"
            CFLAGS="$CFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CFLAGS="$CFLAGS $MFLAGS"
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
            CLIBS="$CLIBS -lstdc++"
#
            CXX='g++-4.3'
            CXXFLAGS='-ansi -D_GNU_SOURCE'
            CXXFLAGS="$CXXFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CXXFLAGS="$CXXFLAGS $MCXXFLAGS $MBAFLAGS -DGLNXA64 -DGCC"
            CXXLIBS="$RPATH $MCXXLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
#
            FC='g95'
            FFLAGS='-fexceptions'
            FFLAGS="$FFLAGS $MFLAGS"
            FLIBS="$RPATH $MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-pthread -shared -Wl,--version-script,$TMW_ROOT/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined"
	    LDFLAGS="$LDFLAGS"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
	    LDEXTENSION='.mexa64'
#
            POSTLINK_CMDS=''
#----------------------------------------------------------------------------
            ;;
        sol64)
#----------------------------------------------------------------------------
            CC='cc -xarch=v9a'
            CFLAGS='-dalign -xlibmieee -D__EXTENSIONS__ -D_POSIX_C_SOURCE=199506L -mt'
            CFLAGS="$CFLAGS $MFLAGS"
            CLIBS="$MLIBS -lm"
            COPTIMFLAGS='-xO3 -xlibmil -DNDEBUG'
            CDEBUGFLAGS='-g'  
#           
            CXX='CC -xarch=v9a -compat=5'
            CCV=`CC -xarch=v9a -V 2>&1`
            version=`expr "$CCV" : '.*\([0-9][0-9]*\)\.'`
            if [ "$version" = "4" ]; then
                    echo "SC5.0 or later C++ compiler is required"
            fi
            CXXFLAGS='-dalign -xlibmieee -D__EXTENSIONS__ -library=stlport4,Crun'
            CXXFLAGS="$CXXFLAGS -D_POSIX_C_SOURCE=199506L -mt"
            CXXFLAGS="$CXXFLAGS $MCXXFLAGS -DSOL64 +d -D_STLP_USE_DECLSPEC=1 -D_STLP_IMPORT_DECLSPEC=__global -D_STLP_CLASS_IMPORT_DECLSPEC=__global"
            CXXLIBS="$MCXXLIBS -library=stlport4,Crun -lm"
            CXXOPTIMFLAGS='-xO3 -xlibmil -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            FC='f90 -xarch=v9a'
            FFLAGS='-dalign -f77=backslash'
            FFLAGS="$FFLAGS $MFLAGS"
            FLIBS="$MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS=''
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'  
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        mac)
#----------------------------------------------------------------------------
echo "Error: Did not imbed 'options.sh' code"; exit 1 #imbed options.sh mac 12
#----------------------------------------------------------------------------
            ;;
        maci)
#----------------------------------------------------------------------------
            CC='gcc-4.0'
            SDKROOT='/Developer/SDKs/MacOSX10.5.sdk'
            MACOSX_DEPLOYMENT_TARGET='10.5'
            ARCHS='i386'
            CFLAGS="-fno-common -no-cpp-precomp -arch $ARCHS -isysroot $SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            CFLAGS="$CFLAGS  -fexceptions"
            CFLAGS="$CFLAGS $MFLAGS"
            CLIBS="$MLIBS"
            COPTIMFLAGS='-O2 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CLIBS="$CLIBS -lstdc++"
            CXX=g++-4.0
            CXXFLAGS="-fno-common -no-cpp-precomp -fexceptions -arch $ARCHS -isysroot $SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            CXXFLAGS="$CXXFLAGS $MCXXFLAGS -DMACI"
            CXXLIBS="$MLIBS -lstdc++"
            CXXOPTIMFLAGS='-O2 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            FC='gfortran'
            FFLAGS='-fexceptions -fbackslash'
            FFLAGS="$FFLAGS $MFLAGS"
            FC_LIBDIR=`$FC -print-file-name=libgfortran.dylib 2>&1 | sed -n '1s/\/*libgfortran\.dylib//p'`
            FC_LIBDIR2=`$FC -print-file-name=libgfortranbegin.a 2>&1 | sed -n '1s/\/*libgfortranbegin\.a//p'`
            FLIBS="$MLIBS -L$FC_LIBDIR -lgfortran -L$FC_LIBDIR2 -lgfortranbegin"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-gdwarf-2'
#
            LD="$CC"
            LDFLAGS="-Wl,-twolevel_namespace -undefined error -arch $ARCHS -Wl,-syslibroot,$SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        maci64)
#----------------------------------------------------------------------------
            CC='gcc-4.0'
            SDKROOT='/Developer/SDKs/MacOSX10.5.sdk'
            MACOSX_DEPLOYMENT_TARGET='10.5'
            ARCHS='x86_64'
            CFLAGS="-fno-common -no-cpp-precomp -arch $ARCHS -isysroot $SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            CFLAGS="$CFLAGS  -fexceptions"
            CFLAGS="$CFLAGS $MFLAGS"
            CLIBS="$MLIBS"
            COPTIMFLAGS='-O2 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CLIBS="$CLIBS -lstdc++"
            CXX=g++-4.0
            CXXFLAGS="-fno-common -no-cpp-precomp -fexceptions -arch $ARCHS -isysroot $SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            CXXFLAGS="$CXXFLAGS $MCXXFLAGS -DMACI64"
            CXXLIBS="$MLIBS -lstdc++"
            CXXOPTIMFLAGS='-O2 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            FC='gfortran'
            FFLAGS='-fexceptions -m64 -fbackslash'
            FFLAGS="$FFLAGS $MFLAGS"
            FC_LIBDIR=`$FC -print-file-name=libgfortran.dylib 2>&1 | sed -n '1s/\/*libgfortran\.dylib//p'`
            FC_LIBDIR2=`$FC -print-file-name=libgfortranbegin.a 2>&1 | sed -n '1s/\/*libgfortranbegin\.a//p'`
            FLIBS="$MLIBS -L$FC_LIBDIR -lgfortran -L$FC_LIBDIR2 -lgfortranbegin"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$CC"
            LDFLAGS="-Wl,-twolevel_namespace -undefined error -arch $ARCHS -Wl,-syslibroot,$SDKROOT -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#           LDEXTENSION="$LDEXTENSION"
#----------------------------------------------------------------------------
#############################################################################
