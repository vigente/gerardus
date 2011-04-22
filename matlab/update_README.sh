#!/bin/bash

# Script to automatically extract the summary of the help headers of
# all the Matlab functions in this toolbox, and generate a README file

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
for DIR in `find . -maxdepth 1 -name "*Toolbox" | sort` 
do

    echo "$DIR" | tr -d './'
    echo '-------------------------------------------------------------'
    echo ''
    
    # loop every function
    for FILE in `find $DIR/*.m | sort`
    do
	echo `basename "$FILE"`
	echo ''
	# get first text block in the header
	# remove the line(s) that declares the function
	# remove the comment characters %
	# keep only the summary of the help header, not the syntax
	grep -m 1 -B 1000  "^$" "$FILE" \
	    | grep '%' \
	    | tr -d '%' \
	    | grep -m 1 -B 100  "^$" \
	    | sed 's/^/\t/'
    done

    echo ''

done
} >> README
