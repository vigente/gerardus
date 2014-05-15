#!/bin/sh
# $Id: WB08_Fig_5.sh,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
#
# This creates Figure 5 which depends on data produced
# by the corresponding Matlab/Octave script.
#
# Wessel, P. and J. M. Becker, 2008, Interpolation using a
#  generalized Green's function for a spherical surface spline
#  in tension, Geophys. J. Int., doi:10.1111/j.1365-246X.2008.03829.x

# Example of gridding a larger spherical data set.  Here we use
# lidar topography measurements from the Moon; see reference
# D. E. Smith, M. T. Zuber, G. A. Neumann and F. G. Lemoine,
# Topography of the Moon from the Clementine lidar. J. Geophys. Res.,
# 102 (1997), pp. 1591–1611.

PS=`basename $0 '.sh'`.ps
PDF=`basename $0 '.sh'`.pdf

# Run WB08_Fig_5.m in Matlab first, then the xyz2grd calls:
p=25
xyz2grd Fig_5a.d -Rg -I15m -Glun2_p0.grd -V
xyz2grd Fig_5b.d -Rg -I15m -Glun2_p${p}.grd -V

makecpt -Crainbow -T-11000/9000/500 -Z > t.cpt
grdgradient lun2_p0.grd -Nt0.5 -M -A45 -Glun2_p0_int.grd
grdgradient lun2_p${p}.grd -Nt1 -M -A45 -Glun2_p${p}_int.grd
grdimage lun2_p${p}.grd -Ilun2_p${p}_int.grd -Ct.cpt -Ei -B30g30 -JH180/6i -P -K -U/-0.5i/-1.75i/"Wessel and Becker: Figure 5" -X0.75i -Y2i > $PS
echo "0 90 14 0 1 LB b)" | pstext -R -J -O -K -N -D-3i/-0.2i >> $PS
psxy -R -J -O -K -Sc0.025i -B30g30 -Gblack lun2.xyz -Y3.25i >> $PS
echo "0 90 14 0 1 LB a)" | pstext -R -J -O -K -N -D-3i/-0.2i >> $PS
psxy -R -J -O /dev/null >> $PS
ps2raster -Tf $PS
ps2raster -Au -Te $PS
open $PDF
rm -f $PS
