#!/bin/bash

# add_filter_template.sh - Script that partly automates adding a new
# filter to itk_imfilter().
#
# If the new filter is part of ITK, it doesn't take user-provided
# input parameters, and returns a single output of the same type as
# the input image, and doesn't require setup steps, then running this
# script will suffice to add the filter to itk_imfilter().
#
# Syntax:
#
#  $ ./add_filter_template.sh XXXImageFilter shortname
#
# Example:
#
#  $ ./add_filter_template.sh AnisotropicDiffusionVesselEnhancementImageFilter advess

# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2012 University of Oxford
# Version: 0.1.3
# $Rev$
# $Date$
#
# University of Oxford means the Chancellor, Masters and Scholars of
# the University of Oxford, having an administrative office at
# Wellington Square, Oxford OX1 2JD, UK. 
#
# This file is part of Gerardus.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. The offer of this
# program under the terms of the License is subject to the License
# being interpreted in accordance with English Law and subject to any
# action against the University of Oxford being under the jurisdiction
# of the English Courts.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

# check syntax
if [ $# -ne 2 ]
    then
    echo "Syntax: `basename $0` XXXImageFilter shortname"
    exit 0
fi  

# mnemonic for input arguments
filter="$1"
shortname="$2"
echo -e "Creating template files to add filter to itk_imfilter():\n    $filter"

# uppercase name of filter for #define declarations
defilter="${filter^^}"

#############################################################
## MexFilter .hpp and .cpp files
#############################################################

pushd MexFilter

# copy template .hpp and .cpp files
cp MexTemplateFilter.hpp Mex${filter}.hpp
cp MexTemplateFilter.cpp Mex${filter}.cpp

# replace template strings by others specific to the new filter
sed -i s/TemplateImageFilter/"$filter"/g \
    Mex${filter}.hpp Mex${filter}.cpp
sed -i s/TEMPLATEIMAGEFILTER/"$defilter"/g \
    Mex${filter}.hpp Mex${filter}.cpp
sed -i -e "s/Copyright.*/Copyright © `date +%Y` University of Oxford/g" \
    Mex${filter}.hpp Mex${filter}.cpp
sed -i 's/Version:.*/Version: 0.1.0/g' \
    Mex${filter}.hpp Mex${filter}.cpp
sed -i 's/$Rev:.*/$Rev$/g' \
    Mex${filter}.hpp Mex${filter}.cpp
sed -i 's/$Date:.*/$Date$/g' \
    Mex${filter}.hpp Mex${filter}.cpp
sed -i -e "s/shortname = \"template\"/shortname = \"${shortname}\"/g" \
    Mex${filter}.hpp Mex${filter}.cpp

echo "    ... Created Mex${filter}.hpp/cpp"

popd

#############################################################
## ItkImFilter.cpp
#############################################################

# #include MexFilter.hpp
line="#include \"Mex${filter}.hpp\""
if [ `grep -c "$line" ItkImFilter.cpp` == 0 ]
then
    sed -i "/#include \"NrrdImage.hpp\"/i $line" \
	ItkImFilter.cpp
fi

# add new filter to enum list of supported filters
line="nMex${filter},"
if [ `grep -c "$line" ItkImFilter.cpp` == 0 ]
then
    sed -i "/enum SupportedFilter {/a \  ${line}" \
	ItkImFilter.cpp
fi

# add macro to instantiate new filter if the corresponding enum is
# selected
line="SELECTFILTER(nMex${filter},\n             Mex${filter})"
if [ `grep -c "SELECTFILTER(nMex${filter}," ItkImFilter.cpp` == 0 ]
then
    sed -i "/#undef SELECTFILTER/i ${line}" \
	ItkImFilter.cpp
fi

# add parsing of filter type to enum
read -r -d '' line <<BLOCK
\
  } else if (ISFILTER(filterName, Mex${filter})) {\n\
\n\
    parseInputTypeToTemplate<nMex${filter}>(nargin, argIn,\n\
							 nargout, argOut);\n
BLOCK
if [ `grep -c "else if (ISFILTER(filterName, Mex${filter}))" ItkImFilter.cpp` == 0 ]
then
    sed -ie "/Insertion point: parseFilterTypeToTemplate/i \  ${line}" \
	ItkImFilter.cpp
fi

echo "    ... Modified ItkImFilter.cpp"

#############################################################
## CmakeLists.txt
#############################################################

# add MexFilter.cpp file to list of source files to build
# itk_imfilter() from
line="Mex${filter}.cpp"
if [ `grep -c "$line" CMakeLists.txt` == 0 ]
then
    sed -i "/  ItkImFilter.cpp/i \  ${line}" \
	CMakeLists.txt
fi

# if the filter is a third-party filter, we need to add its directory
# to the list of includes
pushd ../../cpp/src/third-party > /dev/null
incpath=`find . -name itk${filter}.h`
if [ -n $incpath ]
then
    # remove the file name from the path, keeping only the directory
    incpath=`dirname "$incpath"`
    # remove the trailing ./ from the path
    incpath=${incpath:2}
    # line to add
    line='${GERARDUS_SOURCE_DIR}/cpp/src/third-party/'"${incpath}"
    if [ `grep -c "cpp/src/third-party/$incpath" ../../../matlab/ItkToolbox/CMakeLists.txt` == 0 ]
    then
	sed -ie "/INCLUDE_DIRECTORIES/a \  ${line}" \
	    ../../../matlab/ItkToolbox/CMakeLists.txt
    fi
fi
popd > /dev/null

echo "    ... Modified CMakeLists.txt"
