#!/bin/sh
# $Id: WB08_Fig_2.sh,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
#
# This creates Figure 2 which depends on data produced
# by the corresponding Matlab/Octave script.
#
# Wessel, P. and J. M. Becker, 2008, Interpolation using a
#  generalized Green's function for a spherical surface spline
#  in tension, Geophys. J. Int., doi:10.1111/j.1365-246X.2008.03829.x

PS=`basename $0 '.sh'`.ps
PDF=`basename $0 '.sh'`.pdf

xyz2grd -R0/360/0/90 -I1 -GFig_2_p0.grd Fig_2_p0.d
xyz2grd -R0/360/0/90 -I1 -GFig_2_orig.grd Fig_2_orig.d
pscoast -R0/360/0/90 -JA0/90/6i -P -Glightgray -K -B30 \
	-U/-0.75i/-1.75i/"Wessel and Becker: Figure 2" -Y2i --FRAME_WIDTH=0.025i \
	--PLOT_DEGREE_FORMAT=ddd:mm:ssF --ANNOT_FONT_SIZE=10 > $PS
echo 0 90 | psxy -R -J -O -K -Sx0.1i -W0.5p >> $PS
grdcontour -R Fig_2_p0.grd -J -O -K -Z0.001 -C5 -A10 -Gl195/0/0/90,0/90/295/0 >> $PS
grdcontour -R Fig_2_orig.grd -J -O -K -Z0.001 -C5 -A10 -Wa0.75p,- -Wc0.25p,. -Gl160/0/270/80,270/80/340/0 >> $PS
psxy -R -J -O -K mag_obs_1990.d -Sc0.1 -Gblack >> $PS
psxy -R -J -O -K mag_obs_1990.d -Sc0.025 -Gwhite >> $PS
psxy -R -J -O -K mag_validate_1990.d -Sc0.1 -Gwhite -W0.25p  >> $PS
psxy -R -J -O -K mag_validate_1990.d -Sc0.025 -Gblack >> $PS
echo 104.50 59.92 | psxy -R -J -O -K -Sc0.1 -Gwhite -W0.25p >> $PS
echo 104.50 59.92 | psxy -R -J -O -K -Sx0.1 -W1p >> $PS
echo 104.50 59.92 -42 1.65 | psxy -R -J -O -K -SV0.02i/0.1i/0.08i -Gblack --VECTOR_SHAPE=0.5 >> $PS
psxy -R -J -O /dev/null >> $PS

ps2raster -Tf $PS
ps2raster -Au -Te $PS
open $PDF
rm -f $PS
