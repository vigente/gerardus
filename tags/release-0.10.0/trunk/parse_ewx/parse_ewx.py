#!/usr/bin/env python

# Copyright (c) The Chancellor, Masters and Scholars of the University of
# Oxford, 15th March 2010. All rights reserved.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
A basic parser for EvokeWave (http://www.unisip.hu/evokewave/evokewave.htm)
data files.

It will extract the primary data section to a folder full of text-format files,
suitable for plotting in Gnuplot.  It expects the data to consist of a sequence
of action potentials, with 'junk data' (of amplitude roughly zero) inbetween.
Each action potential will be stored in a separate file (numbered sequentially
starting from zero) with the voltage values listed one per line.

No time data is currently read from the input file, so you have to know what
values to use for the time axis.

The data is stored in a packed binary format of single-precision floating point
values.  The parser uses a simple state machine to track what is an AP and what
is junk.  There are 3 states:
 * No_data - the initial state, and also the state while reading 'junk data'.
   In the latter case we can enter this state from AP_resting iff we see a
   data entry with absolute value less than the --rubbish-threshold.
 * AP_resting - this state is entered whenever we see a data entry with value
   below --v-threshold.
 * AP_triggered - this state is entered from AP_resting if we see a value
   above --v-threshold that doesn't trigger the No_data state.

If real APs are missed, or junk data is included, then try tweaking the two
threshold parameters.
"""

import optparse
import os
import struct
import sys

usage = '%prog [options] <file.ewx> <output_dir>\n' + __doc__
parser = optparse.OptionParser(usage=usage)
parser.add_option('-V', '--v-threshold',
                  type='float', default=-20.0,
                  help="Threshold for identifying the 'at rest' state."
                  " [default: %default]")
parser.add_option('-r', '--rubbish-threshold',
                  type='float', default=1.0,
                  help="Threshold for identifying rubbish data."
                  " [default: %default]")
parser.add_option('-o', '--output-name-format',
                  default='ap%06d.dat',
                  help="Format string for names of output files."
                  " [default: %default]")
parser.add_option('-n', '--stop-after-n-data',
                  type='int', default=0, metavar='N',
                  help="Don't parse the whole file; stop after N data points")
parser.add_option('-a', '--stop-after-n-aps',
                  type='int', default=0, metavar='N',
                  help="Don't parse the whole file; stop after N APs")

options, args = parser.parse_args(sys.argv[1:])
if len(args) != 2:
    parser.error("You need to specify both an input file and output directory.")
fname = args[0]
outdir = args[1]

# Open data file and create output folder
fp = file(fname, 'rb') # Note: must open in binary format!
if os.path.exists(outdir):
    parser.error("Output directory exists.  Please remove it before running.")
os.mkdir(outdir)

# Ignore header lines
for line in fp:
    if line.startswith('Data_'):
        break

def open_output_file(file_num):
    """Create a new output file."""
    return open(os.path.join(outdir, options.output_name_format % file_num),
                'w')

def parse_data(name, fp, ending):
    """The actual parsing function."""
    items_done = 0   # How many data items we've seen
    at_rest = True   # Whether we're in the AP_resting state
    real_data = True # Whether we're not in the No_data state
    file_num = 0     # How many APs we've seen
    outfile = open_output_file(file_num)
    
    print "Parsing", name, "data section"
    data, offset = "", 1
    for line in fp:
        if line.startswith(ending):
            break
        data = data + line
        num_floats = (len(data)-1)//4
        for f in struct.unpack_from('%df' % num_floats, data, offset):
            if at_rest and abs(f) < options.rubbish_threshold:
                # Transit to No_data.
                # Assume we've finished a trace; go to next file
                at_rest = False
                real_data = False
                file_num += 1
                if options.stop_after_n_aps and file_num >= options.stop_after_n_aps:
                    return
                outfile = open_output_file(file_num)
            elif f < options.v_threshold:
                # Transit to AP_resting.
                real_data = True
                at_rest = True
            else:
                # Either transit to AP_triggered, or remain in No_data.
                at_rest = False
            items_done += 1
            if real_data:
                print >>outfile, f
            if options.stop_after_n_data and items_done > options.stop_after_n_data:
                return
        data = data[num_floats*4 + offset:]
        offset = 0


# Primary data section
parse_data('Primary', fp, 'Data_')


## Calced data section
#parse_data('Calced', fp, '[NOTEBOOK]')
