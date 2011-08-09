#!/bin/bash

# add_itk_imfilter_template.sh
#
# Script to automatically create templates and modify ItkImFilter.cpp
# to add a new ITK filter to itk_imfilter(). This script only creates
# the basic skeleton. Specific parameters, methods, etc for the filter
# need to be added by hand.
#
# Syntax:
#
#   $ ./add_itk_imfilter_template.sh file [file]
#
# Example
#
#   $ ./add_itk_imfilter_template.sh itkBinaryErodeImageFilter.h itkBinaryDilateImageFilter.h
#

# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2011 University of Oxford
# Version: 0.2.0
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

# syntax
if [ $# -lt 1 ]
    then
    echo "Syntax: `basename $0` FilterName [FilterName ...]"
    exit 0
fi  

# constants
ITKINCDIR=/usr/local/include/InsightToolkit/
SKIPPEDLOG=add_itk_imfilter_template.skipped
ADDEDLOG=add_itk_imfilter_template.added

# clear log files
rm -f "$SKIPPEDLOG" "$ADDEDLOG"

# auxiliary file to concatenate other files
TMPFILE=`mktemp add_itk_imfilter_template.XXXXXXXXXX`
# log file with the message for the svn commit
LOGFILE=`mktemp add_itk_imfilter_template.XXXXXXXXXX`

for filter in "$@"
do
    echo "$filter" ...
    FILTERPATH=`find "$ITKINCDIR" -name "$filter"`
    if [ -z "$FILTERPATH" ]
    then
	# filter not found, add to the log and skip
	echo -e "\t... not found in ITK (skipping)"
	echo "$filter" >> "$SKIPPEDLOG"
	continue
    else
	# add entry to the found log
	echo -e "\t... found in ITK"
    fi

    # basename of templates without the extension
    MEXFILTERFILE=Mex"${filter##itk}"
    MEXFILTERFILE="${MEXFILTERFILE%.*}"

    # check whether the filter files already exist in Gerardus
    if [ -e ItkToolbox/"$MEXFILTERFILE".hpp -o \
	-e ItkToolbox/"$MEXFILTERFILE".cpp ]
    then
	echo -e "\t... found in Gerardus (skipping)"
	echo "$filter" >> "$SKIPPEDLOG"
	continue
    else
	echo -e "\t... not found in Gerardus"
    fi

    #####################################################################
    ## MexBinaryDilateImageFilter.hpp
    ## MexBinaryDilateImageFilter.cpp
    #####################################################################

    # create new .hpp and .cpp template files for this filter based
    # on the filter template
    svn cp ItkToolbox/MexTemplateFilter.hpp ItkToolbox/"$MEXFILTERFILE".hpp
    svn cp ItkToolbox/MexTemplateFilter.cpp ItkToolbox/"$MEXFILTERFILE".cpp

    echo -e \
	"\t* Add ItkToolbox/"$MEXFILTERFILE".hpp,\n\tItkToolbox/"$MEXFILTERFILE".cpp:" \
	>> $LOGFILE
    echo >> $LOGFILE
    echo -e \
	"\t- Copied from MexTemplateFilter.* to add support for\n\t$filter to itk_imfilter()" \
	>> $LOGFILE
    echo >> $LOGFILE

    # replace generic "template" tags in the files with the specific
    # tags for this filter

    # TemplateImageFilter -> BinaryDilateImageFilter
    typeset +u +l TAG=${MEXFILTERFILE#Mex}
    sed -i s/TemplateImageFilter/$TAG/ ItkToolbox/"$MEXFILTERFILE".?pp
    # TEMPLATEIMAGEFILTER -> BINARYDILATEIMAGEFILTER
    typeset -u TAG=$TAG
    sed -i s/TEMPLATEIMAGEFILTER/$TAG/ ItkToolbox/"$MEXFILTERFILE".?pp
    # std::string>::shortname = "template"; -> std::string>::shortname = "binarydilate";
    typeset -l TAG=${TAG%IMAGEFILTER*}
    sed -i s/'"'template'"'/'"'"$TAG"'"'/ ItkToolbox/"$MEXFILTERFILE".cpp

    #####################################################################
    ## ItkImFilter.cpp
    #####################################################################

    # copy first part of the file to temp file
    sed -n '1,/Gerardus headers/ p' ItkToolbox/ItkImFilter.cpp > $TMPFILE

    # get the block of Gerardus headers
    # remove the first and last lines
    # insert the new #include line at the beginning
    # sort in alphabetical order
    # concatenate to the temp file
    sed -n "/Gerardus headers/,/End Gerardus headers/ p" ItkToolbox/ItkImFilter.cpp \
	| tail -n +2 | head -n -1 \
	| sed "1 i #include \"$MEXFILTERFILE.hpp\"" \
	| sort >> $TMPFILE
    
    # append code until next insertion point
    # remove insertion point line
    sed -n '/End Gerardus headers/,/Insertion point: parseFilterTypeAndRun/ p' \
	ItkToolbox/ItkImFilter.cpp \
	| head -n -1 >> $TMPFILE

    # append lines with the "else if" block for the new filter
    echo '  } else if (ISFILTER(filterName, '"$MEXFILTERFILE"')) {' \
	>> $TMPFILE
    echo >> $TMPFILE
    echo '    #error Input arguments cannot be determined automatically by' `basename $0` >> $TMPFILE
    echo '    filter = new '"$MEXFILTERFILE"'<InVoxelType,' >> $TMPFILE
    echo '                    OutVoxelType>(nrrd, nargout, argOut, nargin, argIn);' \
	>> $TMPFILE
    echo >> $TMPFILE

    # append code until next insertion point
    # remove insertion point line
    sed -n '/Insertion point: parseFilterTypeAndRun/,/Insertion point: parseOutputTypeToTemplate/ p' \
	ItkToolbox/ItkImFilter.cpp \
	| head -n -1 >> $TMPFILE

    # append lines with the "else if" block for the new filter
    echo '  } else if (ISFILTER(filter, '"$MEXFILTERFILE"')) {' \
	>> $TMPFILE
    echo >> $TMPFILE
    echo '    #error outVoxelType cannot be determined automatically by' `basename $0` >> $TMPFILE
    echo '    outVoxelType = SAME;' >> $TMPFILE
    echo >> $TMPFILE

    # append rest of the code
    sed -n '/Insertion point: parseOutputTypeToTemplate/,// p' \
	ItkToolbox/ItkImFilter.cpp >> $TMPFILE

    # update ItkImFilter with the new version
    cp $TMPFILE ItkToolbox/ItkImFilter.cpp

    # keep a list of added filters
    echo -e "\t\t$filter" >> "$ADDEDLOG"

done

##################################################
# ChangeLog
##################################################

# today's entry header
DATE=`date +%Y-%m-%d`
echo $DATE'  Ramon Casero  <ramon.casero@cs.ox.ac.uk>' > $TMPFILE
echo >> $TMPFILE

# ItkImFilter change entry
echo -e "\t* ItkToolbox/ItkImFilter.cpp (v0.0.0):" >> $TMPFILE
echo >> $TMPFILE
echo -e "\t- Added support for filters:" >> $TMPFILE
echo >> $TMPFILE
cat "$ADDEDLOG" >> $TMPFILE
echo  >> $TMPFILE

# New filter .cpp .hpp files entry
cat $LOGFILE >> $TMPFILE

# if there's an entry in the ChangeLog for today, we don't
# duplicate the line with the date
TODAYHASANENTRY=`head -n 1 ChangeLog | grep -m 1 -C 0 $DATE`
if [ -z "$TODAYHASANENTRY" ]
then
    cat ChangeLog >> $TMPFILE
else
    tail -n +3 ChangeLog >> $TMPFILE
fi
mv $TMPFILE ChangeLog

# delete temporal files
rm -f $TMPFILE $LOGFILE
