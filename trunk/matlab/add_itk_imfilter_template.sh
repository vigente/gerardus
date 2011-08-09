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
# Version: 0.1.0
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

    # create new .hpp and .cpp template files for this filter based
    # on the filter template
    svn cp ItkToolbox/MexTemplateFilter.hpp ItkToolbox/"$MEXFILTERFILE".hpp
    svn cp ItkToolbox/MexTemplateFilter.cpp ItkToolbox/"$MEXFILTERFILE".cpp

    echo -e \
	"\tAdd ItkToolbox/"$MEXFILTERFILE".hpp,\n\tItkToolbox/"$MEXFILTERFILE".cpp:" \
	>> $LOGFILE
    echo >> $LOGFILE
    echo -e \
	"\t- Copied from MexTemplateFilter.* to add support for\n\t$filter to itk_imfilter()" \
	>> $LOGFILE
    echo >> $LOGFILE

    # filter has been added
    echo ItkToolbox/"$MEXFILTERFILE".hpp >> "$ADDEDLOG"
    echo ItkToolbox/"$MEXFILTERFILE".cpp >> "$ADDEDLOG"
done

# update the ChangeLog
DATE=`date +%Y-%m-%d`
echo $DATE'  Ramon Casero  <ramon.casero@cs.ox.ac.uk>' > $TMPFILE
echo >> $TMPFILE
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

cat "$ADDEDLOG" | xargs svn ci -F "$LOGFILE" ChangeLog 

# delete temporal files
rm -f $TMPFILE $LOGFILE
