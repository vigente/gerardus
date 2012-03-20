#!/bin/bash

# images2avi.sh - Create an AVI video from a stack of images
#
# Usage:
#
#    ./images2avi.sh file1.jpg [file2.jpg...]
#    ./images2avi.sh file1.tif [file2.tif...]
#    ./images2avi.sh file1.png [file2.png...]
#
# The output in the current directory will be:
#
#    * A video file called out.avi.
#
#    * A temp directory with a name with the format images2avi.XXXXXXXX
#      that contains the resized JPEG files from which the video is
#      created.
#
# In principle, any supported image type by ImageMagick can be
# used. File formats can be mixed too, as all formats will be
# converted to JPEG before creating the video.
#
# For the moment, this script does not accept input arguments. If you
# want to change the internal parameters:
#
#    * frame rate=25
#    * name of output file=out.avi
#    * horizontal number of pixels=640).
#
# you need to edit this script.

# Copyright © 2009 University of Oxford
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# PARAMETERS
# size of the horizontal side of each frame in the final video
outsize=640
outfile=out.avi
frameps=25

# create temporal directory for intermediate files
tmpdir=`mktemp -d images2avi.XXXXXXXX` || exit 1

# loop through each file given in the command line as inputs
for file in "$@"
do
    
    # file to process
    echo "$file"

    # get basename without the extension
    filename=`basename "$file"`
    filename="${filename%.*}"

    # get image dimensions
    width=`identify -format "%w" "$file"`
    height=`identify -format "%h" "$file"`

    # convert file to JPEG, in the correct size for the video, 
    # respecting the aspect ratio, and copy to the temp directory
    convert -size "$width"x"$height" -resize "${outsize}x>" \
	    "$file" "$tmpdir/$filename.jpg"


done

# create the video
mencoder "mf://$tmpdir/*.jpg" -mf fps="$frameps" -o "$outfile"  \
    -ovc lavc -lavcopts vcodec=mpeg4

# delete the temporary directory
echo
echo
echo "Delete temp directory with output frames? (y/N)"
read delete
if [ "$delete" == "y" ]
then
    rm -rf "$tmpdir"
fi
