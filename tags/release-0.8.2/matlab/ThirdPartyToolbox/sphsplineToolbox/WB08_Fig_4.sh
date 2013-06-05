#!/bin/sh
# $Id: WB08_Fig_4.sh,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
#
# This creates Figure 4 which depends on data produced
# by the corresponding Matlab/Octave script.
#
# Wessel, P. and J. M. Becker, 2008, Interpolation using a
#  generalized Green's function for a spherical surface spline
#  in tension, Geophys. J. Int., doi:10.1111/j.1365-246X.2008.03829.x

PS=`basename $0 '.sh'`.ps
PDF=`basename $0 '.sh'`.pdf
# This example uses 370 radio occultation data for Mars to grid the topography.
# Data and information from
# Smith, D. E., and M. T. Zuber (1996), The shape of Mars and the topographic
# signature of the hemispheric dichotomy, Science, 271, 184–187.

# Make Mars ellipsoid given their three best-fitting axes:
a=3399.472
b=3394.329
c=3376.502
grdmath -Rg -I1 X COSD $a DIV DUP MUL X SIND $b DIV DUP MUL ADD Y COSD DUP MUL MUL Y SIND $c DIV DUP MUL ADD SQRT INV = ellipsoid.grd

# Run WB08_Fig_4.m in Matlab first, then the xyz2grd calls:
xyz2grd Fig_4a.d -Rg -I1 -Gmars.grd -V
xyz2grd Fig_4b.d -Rg -I1 -Gmars2.grd -V
# Scale to km and remove ellipsoid
grdmath mars.grd 1000 DIV ellipsoid.grd SUB = mars.grd
grdmath mars2.grd 1000 DIV ellipsoid.grd SUB = mars2.grd
makecpt -Crainbow -T-7/15/1 -Z > mars.cpt
grdgradient mars2.grd -M -Ne0.75 -A45 -Gmars2_i.grd
grdimage mars2.grd -Imars2_i.grd -Cmars.cpt -B30g30Wsne -JH0/6i -P -K -X0.75i -Y2i \
	-U/-0.5i/-1.75i/"Wessel and Becker: Figure 4" --ANNOT_FONT_SIZE=12 > $PS
grdcontour mars2.grd -J -O -K -C1 -A5 -Glz+/z- >> $PS
psxy -Rg -J -O -K -Sc0.045i -Gblack mars370.in  >> $PS
echo "0 90 14 0 1 LB b)" | pstext -R -J -O -K -N -D-3i/-0.2i >> $PS
grdgradient mars.grd -M -Ne0.75 -A45 -Gmars_i.grd
grdimage mars.grd -Imars_i.grd -Cmars.cpt -B30g30Wsne -J -O -K -Y3.6i  --ANNOT_FONT_SIZE=12 >> $PS
grdcontour mars.grd -J -O -K -C1 -A5 -Glz+/z- >> $PS
psxy -Rg -J -O -K -Sc0.045i -Gblack mars370.in  >> $PS
psscale -Cmars.cpt -O -K -D3i/-0.1i/5i/0.1ih -I --ANNOT_FONT_SIZE=12 -B2f1/:km: >> $PS
echo "0 90 14 0 1 LB a)" | pstext -R -J -O -N -D-3i/-0.2i >> $PS
ps2raster -Tf $PS
ps2raster -Au -Te $PS
open $PDF
rm -f $PS mars2_i.grd mars_i.grd mars.cpt ellipsoid.grd
