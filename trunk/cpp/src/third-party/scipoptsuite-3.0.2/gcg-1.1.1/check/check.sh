#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       *
#*                         Zuse Institute Berlin (ZIB)                       *
#*                                                                           *
#* This program is free software; you can redistribute it and/or             *
#* modify it under the terms of the GNU Lesser General Public License        *
#* as published by the Free Software Foundation; either version 3            *
#* of the License, or (at your option) any later version.                    *
#*                                                                           *
#* This program is distributed in the hope that it will be useful,           *
#* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
#* GNU Lesser General Public License for more details.                       *
#*                                                                           *
#* You should have received a copy of the GNU Lesser General Public License  *
#* along with this program; if not, write to the Free Software               *
#* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
# @author Martin Bergner
# @author Gerald Gamrath
# @author Christian Puchert

TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
THREADS=$8
FEASTOL=$9
DISPFREQ=${10}
CONTINUE=${11}
LOCK=${12}
VERSION=${13}
LPS=${14}
VALGRIND=${15}
MODE=${16}

SETDIR=../settings

if test ! -e results
then
    mkdir results
fi
if test ! -e locks
then
    mkdir locks
fi

LOCKFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.lock
RUNFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.run.$BINID
DONEFILE=locks/$TSTNAME.$SETNAME.$VERSION.$LPS.done

OUTFILE=results/check.$TSTNAME.$BINID.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINID.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINID.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINID.$SETNAME.tex
TMPFILE=results/check.$TSTNAME.$BINID.$SETNAME.tmp
SETFILE=results/check.$TSTNAME.$BINID.$SETNAME.set

SETTINGS=$SETDIR/$SETNAME.set

if test "$LOCK" = "true"
then
    if test -e $DONEFILE
    then
        echo skipping test due to existing done file $DONEFILE
        exit
    fi
    if test -e $LOCKFILE
    then
        if test -e $RUNFILE
        then
            echo continuing aborted run with run file $RUNFILE
        else
            echo skipping test due to existing lock file $LOCKFILE
            exit
        fi
    fi
    date > $LOCKFILE
    date > $RUNFILE
fi

if test ! -e $OUTFILE
then
    CONTINUE=false
fi

if test "$CONTINUE" = "true"
then
    MVORCP=cp
else
    MVORCP=mv
fi

DATEINT=`date +"%s"`
if test -e $OUTFILE
then
    $MVORCP $OUTFILE $OUTFILE.old-$DATEINT
fi
if test -e $ERRFILE
then
    $MVORCP $ERRFILE $ERRFILE.old-$DATEINT
fi

if test "$CONTINUE" = "true"
then
    LASTPROB=`awk -f getlastprob.awk $OUTFILE`
    echo Continuing benchmark. Last solved instance: $LASTPROB
    echo "" >> $OUTFILE
    echo "----- Continuing from here. Last solved: $LASTPROB -----" >> $OUTFILE
    echo "" >> $OUTFILE
else
    LASTPROB=""
fi

uname -a >>$OUTFILE
uname -a >>$ERRFILE
date >>$OUTFILE
date >>$ERRFILE

# we add 10% to the hard time limit and additional 10 seconds in case of small time limits
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 10\` + \`expr $TIMELIMIT / 10\``

# we add 10% to the hard memory limit and additional 1000mb to the hard memory limit
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 1000\` + \`expr $MEMLIMIT / 10\``
HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024`

echo "hard time limit: $HARDTIMELIMIT s" >>$OUTFILE
echo "hard mem limit: $HARDMEMLIMIT k" >>$OUTFILE

VALGRINDCMD=
if test "$VALGRIND" = "true"
then
   VALGRINDCMD="valgrind --log-fd=1 --leak-check=full"
fi

for i in `cat testset/$TSTNAME.test` DONE
do
    if test "$i" = "DONE"
    then
        date > $DONEFILE
        break
    fi

    if test "$LASTPROB" = ""
    then
        DIR=`dirname $i`
        NAME=`basename $i .gz`
        NAME=`basename $NAME .mps`
        NAME=`basename $NAME .lp`
        BLKFILE=$DIR/$NAME.blk
        DECFILE=$DIR/$NAME.dec
        LASTPROB=""
        if test -f $i
        then
            echo @01 $i ===========
	    NAME=`basename $i`
	    base=${NAME%%.*}
#	    echo $base
            echo @01 $i ===========                >> $ERRFILE
            echo > $TMPFILE
            if test "$SETNAME" != "default"
            then
                echo set load $SETTINGS            >>  $TMPFILE
            fi
            if test "$FEASTOL" != "default"
            then
                echo set numerics feastol $FEASTOL >> $TMPFILE
            fi
            echo set limits time $TIMELIMIT        >> $TMPFILE
            echo set limits nodes $NODELIMIT       >> $TMPFILE
            echo set limits memory $MEMLIMIT       >> $TMPFILE
            echo set lp advanced threads $THREADS  >> $TMPFILE
            echo set timing clocktype 1            >> $TMPFILE
            echo set display verblevel 4           >> $TMPFILE
            echo set display freq $DISPFREQ        >> $TMPFILE
            echo set memory savefac 1.0            >> $TMPFILE # avoid switching to dfs - better abort with memory error
            if test "$LPS" = "none"
            then
                echo set lp solvefreq -1           >> $TMPFILE # avoid solving LPs in case of LPS=none
            fi
            echo set save $SETFILE                 >> $TMPFILE
            echo read $i                           >> $TMPFILE

	    if test $MODE = "detect"
            then
		echo write prob images\/$base.gp  >> $TMPFILE
		echo presolve                      >> $TMPFILE
		echo detect                        >> $TMPFILE
		echo write prob images\/$base-dec.gp  >> $TMPFILE
		echo write prob decs\/$base.dec    >> $TMPFILE
	    elif test $MODE = "detectall"
            then
		echo detect                        >> $TMPFILE
		echo write all ref                 >> $TMPFILE
	    else
		if test $MODE = "readdec"
		then
		    if test -f $DECFILE
		    then
			BLKFILE=$DECFILE
		    fi
		    if test -f $BLKFILE
		    then
			presol=`grep -A1 PRESOLVE $BLKFILE`
		    # if we find a presolving file
			if test $? = 0
			then
                        # look if its in there
			    if grep -xq 1 - <<EOF
$presol
EOF
			    then
				echo presolve          >> $TMPFILE
			    fi
			fi
			echo read $BLKFILE             >> $TMPFILE
		    fi
		fi
		echo optimize                      >> $TMPFILE
		echo display statistics            >> $TMPFILE
#		echo display additionalstatistics  >> $TMPFILE
#               echo display solution                  >> $TMPFILE
		echo checksol                      >> $TMPFILE
	    fi
            echo quit                              >> $TMPFILE
            echo -----------------------------
            date
            date >>$ERRFILE
            echo -----------------------------
            date +"@03 %s"
            bash -c " ulimit -t $HARDTIMELIMIT s; ulimit -v $HARDMEMLIMIT k; ulimit -f 200000; $VALGRINDCMD ../$BINNAME < $TMPFILE" 2>>$ERRFILE
            date +"@04 %s"
            echo -----------------------------
            date
            date >>$ERRFILE
            echo -----------------------------
            echo
            echo =ready=
        else
            echo @02 FILE NOT FOUND: $i ===========
            echo @02 FILE NOT FOUND: $i =========== >>$ERRFILE
        fi
    else
        echo skipping $i
        if test "$LASTPROB" = "$i"
        then
            LASTPROB=""
        fi
    fi
done | tee -a $OUTFILE

rm -f $TMPFILE
rm -f cipreadparsetest.cip

date >>$OUTFILE
date >>$ERRFILE

if test -e $DONEFILE
then
    ./evalcheck.sh $OUTFILE

    if test "$LOCK" = "true"
    then
        rm -f $RUNFILE
    fi
fi
