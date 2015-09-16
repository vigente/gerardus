# imresize-meta.sh
#
# Resize one or more 2D images, taking care of also changing the pixel
# spacing accordingly so that the real world print size remains the
# same.
#
# Original images are not overwritten. Instead, the resized copy is
# saved to a directory chosen by the user.
#
# Syntax:
#
# ./imresize-meta.sh [-o outdir] [-s scale] file1 [file2...]
#
# Options:
#
#   -o outdir: (def /tmp) Output directory to save resized images to
#
#   -s scale:  (def 5) Resizing factor in %. E.g. -s 10 means resize
#              images to 10% of original size.
#
#   file1 [file2...] List of image files to be resized (e.g. PNG, TIF...).

# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2015 University of Oxford
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

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# initialize our own variables:
outdir="/tmp"
verbose=0
scale=5

# read input options
while getopts "vo:s:" opt; do
    case "$opt" in
    v)  verbose=1
        ;;
    s)  scale=$OPTARG
	;;
    o)  outdir=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if [ ${verbose} -eq 1 ]
then
    echo "verbose=$verbose, outdir='$outdir', Leftovers: $@"
fi

# loop files
for file in "$@"
do
    # output file's full path
    outfile="$outdir/${file##*/}"

    if [ ${verbose} -eq 1 ]
    then
	echo "INPUT: $file"
	echo "OUTPUT: $outfile"
    fi

    # get resolution of original image
    x0=`identify -ping -format '%[resolution.x]' "$file"`
    y0=`identify -ping -format '%[resolution.y]' "$file"`

    # compute spacing of scaled image
    x=`bc -l <<< "${x0}/100*${scale}"`
    y=`bc -l <<< "${y0}/100*${scale}"`

    # correct spacing
    convert -resample "${x}x${y}" "$file" "$outfile"

done
