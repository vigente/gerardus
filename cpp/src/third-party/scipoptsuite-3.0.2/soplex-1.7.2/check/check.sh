#!/bin/sh
# $Id: check.sh,v 1.29 2008/09/03 10:47:48 bzfpfets Exp $
# Parameters
# $1 Name of the test, e.g. netlib (needs netlib.test, netlib.solu)
# $2 Path/Name of the binary, e.g. ../bin/soplex.linux.x86.gnu.opt
# $3 Algorithms to test (1...6), e.g. "1 2 3 4"
# $4 Limits, e.g. -l10000 as time limit.
BINNAME=`basename $2`
TSTNAME=`basename $1 .test`
OUTFILE=check.$TSTNAME.$BINNAME.out
ERRFILE=check.$TSTNAME.$BINNAME.err
RESFILE=check.$TSTNAME.$BINNAME.res
date >$OUTFILE
date >$ERRFILE

# Avoid problems with foreign locales (two separate commands for SunOS)
LANG=C
export LANG

# Determine awk program to use.
AWK=awk
OSTYPE=`uname -s | tr '[:upper:]' '[:lower:]' | sed -e s/cygwin.*/cygwin/ -e s/irix../irix/`

case $OSTYPE in
    osf1)  AWK=gawk ;;
    sunos)  AWK=gawk ;;
    aix)  AWK=gawk ;;
esac

#
for i in `cat $1`
do
    echo @01 $i ===========
    echo @01 $i =========== >>$ERRFILE
    for k in $3
    do
        case $k in
	1)  echo =type= LC
	    opt="" ;;
	2)  echo =type= EC
	    opt="-e" ;;
	3)  echo =type= LR
	    opt="-r" ;;
	4)  echo =type= ER
            opt="-e -r" ;;
	5)  echo =type= LCi
	    opt="-i" ;;
	6)  echo =type= ECi
            opt="-e -i" ;;
	7)  echo =type= LCd
	    opt="-p2" ;;
	8)  echo =type= ECd
	    opt="-e -p2" ;;
	9)  echo =type= LCh
	    opt="-t1" ;;
	10) echo =type= ECh
	    opt="-e -t1" ;;
	11) echo =type= LCm
	    opt="-p1" ;;
	12) echo =type= ECm
	    opt="-e -p1" ;;
	#
	# These settings are used for coverage testing and are therefore unusual.
	# They are chosen such that in combination with the above settings every
        # pricer and ratio tester is exercised at least once.
        #
	13) echo =type= CV1
	    opt="-r -i -p0 -t0 -c1" ;;
	14) echo =type= CV2
	    opt="-e -r -i -p3 -t1 -c2 -s3" ;;
	15) echo =type= CV3
	    opt="-r -i -t0 -c3 -s4" ;;
	#
	# Here we test the new bound flipping ratio test in all four combinations of basis representation and algorithm type.
        #
	16) echo =type= LCb
	    opt="-t3" ;;
        17) echo =type= ECb
            opt="-e -t3" ;;
        18) echo =type= LRb
            opt="-r -t3" ;;
        19) echo =type= ERb
            opt="-e -r -t3" ;;
	20) echo =type= R0
	    opt="-d1e-25 -R1e-9" ;;
	21) echo =type= R1
	    opt="-d1e-25 -R1e-9 -s0 -g0" ;;
	22) echo =type= R2
	    opt="-d1e-25 -r" ;;
	23) echo =type= R3
	    opt="-d1e-250 -R1e-6 -e" ;;
        esac
        $2 $opt -C -q $4 $i 2>>$ERRFILE
        echo =ready=
    done
done | tee -a $OUTFILE
date >>$OUTFILE
date >>$ERRFILE
$AWK -f check.awk $TSTNAME.solu $OUTFILE | tee $RESFILE

