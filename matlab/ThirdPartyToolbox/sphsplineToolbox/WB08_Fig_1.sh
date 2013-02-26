#!/bin/sh
# $Id: WB08_Fig_1.sh,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
#
# This creates Figure 1 which depends on data produced
# by the corresponding Matlab/Octave script.
#
# Wessel, P. and J. M. Becker, 2008, Interpolation using a
#  generalized Green's function for a spherical surface spline
#  in tension, Geophys. J. Int., doi:10.1111/j.1365-246X.2008.03829.x

PS=`basename $0 '.sh'`.ps
PDF=`basename $0 '.sh'`.pdf

# Gradients
psxy Fig_1_gp{1,5,25,100}.d -W0.5p -R0/180/-1/0 -JX4i/3i -P -X2.5i -Y3i --CHAR_ENCODING=ISOLatin1+ \
	-B30f10:,-\\232:/0.2f0.1:"|@~\\321@~g(@~g@~)|":WSne -U/-2.25i/-2.75i/"Wessel and Becker: Figure 1" -K > $PS
psxy Fig_1_gp0.d -W1.5p -R -J -O -K >> $PS
pstext -R -J -O -K << EOF >> $PS
123	-0.65	9	0	0	CB 0
116	-0.61	9	0	0	CB 1
75	-0.38	9	0	0	CB 5
38	-0.2	9	0	0	CB 25
23	-0.13	9	0	0	CB 100
EOF
# Surface
psxy Fig_1_p{1,5,25,100}.d -W0.5p -R0/180/0/1 -J -O -K -Y3.25i --CHAR_ENCODING=ISOLatin1+ \
	-B30f10:,-\\232:/0.2f0.1:"g(@~g@~)":WsNe >> $PS
psxy Fig_1_p0.d -W1.5p -R -J -O -K >> $PS
pstext -R -J -O << EOF >> $PS
75	0.5	9	0	0	CB 0
64	0.44	9	0	0	CB 1
45	0.38	9	0	0	CB 5
30	0.3	9	0	0	CB 25
20	0.23	9	0	0	CB 100
EOF
ps2raster -Tf $PS
ps2raster -Au -Te $PS
open $PDF
rm -f $PS
