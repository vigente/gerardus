#!/bin/bash

# Script to automatically extract the summary of the help headers of
# all the Matlab functions in this toolbox, and generate a README file

# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2011 University of Oxford
# Version: 0.1.4
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

# header

{
echo 'This file is part of the Matlab Toolboxes of project Gerardus.'
echo 
echo '============================================================='
echo 'Toolboxes in Gerardus:'
echo '============================================================='
echo ''
echo '* CardiacToolbox'
echo ''
echo '	Functions specific to cardiac image processing'
echo ''
echo '* FileFormatToolbox'
echo ''
echo '	Functions to create image files or convert image files from'
echo '	one format to another.'
echo ''
echo '* FiltersToolbox'
echo ''
echo '	Filters to enhance or transform images in general, and SCI'
echo '	NRRD data volumes in particular.'
echo ''
echo '* ItkToolbox'
echo ''
echo '	ITK functions running as MEX files in Matlab.'
echo ''
echo '* PointsToolbox'
echo ''
echo '	Functions to operate with sets of points.'
echo ''
echo '* ThirdPartyToolbox'
echo ''
echo '	Derivative works or third party functions that cannot be'
echo '	covered by the GPL used elsewhere in Gerardus, or code with an'
echo '	uncertain licence status.'
echo ''
echo ''
} > README

# loop every toolbox
{
for DIR in `find . -maxdepth 2 -name "*Toolbox" | grep -v '/bin' | sort` 
do

    echo "$DIR" | sed 's/^.\///'
    echo '-------------------------------------------------------------'
    echo ''
    
    # loop every function
    for FILE in `find $DIR/*.m | sort`
    do

	case "${DIR}" in
	    
	    # the Spharm Toolbox has a different convention for the
	    # help headers. It's easier to add a special case to this
	    # script than modifiying all the help headers in Spahrm
	    ./ThirdPartyToolbox/SpharmToolbox)

		echo `basename "$FILE"`
		echo ''
		# get the line that says "% goal:" or "% Goal"
		# remove the comment characters %
		# add a tabulation before each line
		# do a line return
		grep -m 1 '% [Gg]oal:' "$FILE" \
		    | tr -d '%' \
		    | sed 's/^/\t/'
		echo ''
		;;
	    
	    
	    # the iso2mesh Toolbox has a different convention for the
	    # help headers. It's easier to add a special case to this
	    # script than modifiying all the help headers in iso2mesh
	    ./ThirdPartyToolbox/iso2meshToolbox)

		echo `basename "$FILE"`
		echo ''
        	# get first text block in the header
                # remove the line(s) that declares the function
        	# remove the comment characters %
        	# keep only the summary of the help header, not the syntax
		# add a tabulation before each line
		grep -m 1 -B 1000  "^$" "$FILE" \
		    | grep '%' \
		    | tr -d '%' \
		    | tail -n +4 \
		    | grep -m 1 -B 100  "^$" \
		    | sed 's/^/\t/'
		;;
	    
	    *)
		echo `basename "$FILE"`
		echo ''
        	# get first text block in the header
                # remove the line(s) that declares the function
        	# remove the comment characters %
        	# keep only the summary of the help header, not the syntax
		# add a tabulation before each line
		grep -m 1 -B 1000  "^$" "$FILE" \
		    | grep '%' \
		    | tr -d '%' \
		    | grep -m 1 -B 100  "^$" \
		    | sed 's/^/\t/'
		;;
	esac
    done
    
    echo ''

done
} >> README
